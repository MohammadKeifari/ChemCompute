from dataclasses import dataclass
import math
from typing import Optional, Union

import numpy as np

METHOD_DEFAULTS = {
    "bgd": {"max_iter": 5000, "learning_rate": 0.1, "tol": 1e-8},
    "sgd": {"max_iter": 5000, "learning_rate": 0.1, "tol": 1e-8},
    "newton": {"max_iter": 200, "learning_rate": 1.0, "tol": 1e-10},
}

VALID_METHODS = tuple(METHOD_DEFAULTS.keys())
INFINITE_K_VALUE = 1e300
INFINITE_LNK = math.log(INFINITE_K_VALUE)


def _huber_value(x: np.ndarray, delta: float) -> np.ndarray:
    ax = np.abs(x)
    return np.where(ax <= delta, 0.5 * x ** 2, delta * (ax - 0.5 * delta))


def _huber_derivative(x: np.ndarray, delta: float) -> np.ndarray:
    ax = np.abs(x)
    return np.where(ax <= delta, x, delta * np.sign(x))


@dataclass
class EquilibriumLoss:
    name: str
    delta: float = 1.0

    def residual(self, lnQ: np.ndarray, lnK: np.ndarray) -> np.ndarray:
        log_diff = lnQ - lnK
        if self.name == "log_quotient":
            return log_diff
        if self.name == "quotient_error":
            return np.expm1(log_diff)
        if self.name == "log_huber":
            return log_diff
        raise ValueError(f"Unknown loss: {self.name}")

    def objective(self, residual: np.ndarray) -> float:
        if self.name == "log_huber":
            return float(np.sum(_huber_value(residual, self.delta)))
        return float(0.5 * np.dot(residual, residual))

    def grad_weights(self, residual: np.ndarray) -> np.ndarray:
        if self.name == "log_huber":
            return _huber_derivative(residual, self.delta)
        return residual

    def jacobian_scale(self, residual: np.ndarray) -> np.ndarray:
        if self.name == "quotient_error":
            return 1.0 + residual
        return np.ones_like(residual)


LOSS_REGISTRY = {
    "log_quotient": EquilibriumLoss("log_quotient"),
    "quotient_error": EquilibriumLoss("quotient_error"),
    "log_huber": EquilibriumLoss("log_huber"),
}


@dataclass
class EquilibriumResult:
    """Equilibrium solve output and per-reaction diagnostics."""

    concentrations: list[float]
    compounds: list[str]
    reaction_extents: list[float]
    """Stoichiometric extent x per reaction (moles advanced); informational only."""
    reaction_quotient_error: list[float]
    """Per-reaction |Q/K - 1| (0 for infinite-K reactions). A value of 1.0 means Q/K ≈ 0 or 2."""
    max_reaction_quotient_error: float
    reaction_quotient_ratio: list[float]
    converged: bool
    stop_reason: str
    iterations: int
    criterion_met: bool
    criterion_type: str
    criterion_value: float
    criterion_limit: float
    electrode_Eh: Optional[float] = None
    """Imposed or solved electrode potential (V vs SHE), when redox coupling is active."""

    @property
    def q_over_k(self) -> list[float]:
        """Per-reaction Q/K ratio at the solution (1.0 means equilibrium)."""
        return self.reaction_quotient_ratio

    @property
    def concentrations_dict(self) -> dict[str, float]:
        """Map compound formula labels to equilibrium concentrations."""
        return dict(zip(self.compounds, self.concentrations))


@dataclass
class EquilibriumContext:
    N: np.ndarray
    S: np.ndarray
    A: np.ndarray
    c0: np.ndarray
    lnK: np.ndarray
    R: int
    C: int
    min_concentration: float
    infinite_k_mask: np.ndarray
    excess_source_mask: np.ndarray
    env: object = None
    ln_gamma: np.ndarray = None
    redox_indices: np.ndarray = None
    couple_eh: bool = False


def _update_activity(ctx: EquilibriumContext, c: np.ndarray) -> None:
    """Refresh ln(gamma) from current concentrations when activity model is active."""
    if ctx.env is None or getattr(ctx.env, "activity_model", None) is None:
        ctx.ln_gamma = np.zeros(ctx.C, dtype=float)
        return
    gammas = ctx.env.activity_model.gamma_array(ctx.env, np.maximum(c, 0.0), ctx.env.T)
    ctx.ln_gamma = np.log(np.maximum(gammas, 1e-30))


def _infinite_k_extent(ctx: EquilibriumContext, reaction_index: int, x: np.ndarray) -> float:
    """Additional forward extent for an irreversible reaction from the current state."""
    c = np.maximum(ctx.c0 + ctx.S @ x, 0.0)
    lower = -np.inf
    upper = np.inf
    for j in range(ctx.C):
        s = ctx.S[j, reaction_index]
        if abs(s) < 1e-15:
            continue
        bound = -c[j] / s
        if s > 0:
            lower = max(lower, bound)
        else:
            upper = min(upper, bound)
    if lower == -np.inf or lower > upper:
        return x[reaction_index]
    return x[reaction_index] + lower


def _apply_infinite_k_extents(ctx: EquilibriumContext, x: np.ndarray) -> np.ndarray:
    if not ctx.infinite_k_mask.any():
        return x
    for _ in range(ctx.R):
        changed = False
        for r in range(ctx.R):
            if not ctx.infinite_k_mask[r]:
                continue
            new_value = _infinite_k_extent(ctx, r, x)
            if not np.isclose(new_value, x[r], rtol=0.0, atol=1e-15):
                x[r] = new_value
                changed = True
        if not changed:
            break
    return x


def _initialize_extents(ctx: EquilibriumContext) -> np.ndarray:
    x = np.zeros(ctx.R, dtype=float)
    return _apply_infinite_k_extents(ctx, x)


def _safe_concentrations(c: np.ndarray, min_concentration: float) -> np.ndarray:
    """Use the true positive concentration for logs; floor only non-positive values."""
    return np.where(c > 0, c, min_concentration)


def _build_context(env, min_concentration: float) -> EquilibriumContext:
    N = env.stoichiometric_coefficient_array
    S = N.T.copy()
    R, C = N.shape

    c0 = np.array(env.concentrations, dtype=float)
    A = -N.astype(float)

    K_vec = np.empty(R, dtype=float)
    for i, rxn in enumerate(env.reactions):
        if getattr(rxn, "infinite_K", False):
            K_vec[i] = INFINITE_K_VALUE
        else:
            K_vec[i] = max(rxn.K, 1e-300)
    lnK = np.log(K_vec)
    infinite_k_mask = np.array([getattr(rxn, "infinite_K", False) for rxn in env.reactions], dtype=bool)
    excess_source_mask = np.zeros(R, dtype=bool)
    for i, rxn in enumerate(env.reactions):
        for species in rxn.reactants + rxn.products:
            compound = species["compound"]
            if species.get("excess", False) or getattr(compound, "excess", False):
                excess_source_mask[i] = True
                break

    for j, compound in enumerate(env.compounds):
        ph = compound.phase(env.T)
        if ph in ("s", "l"):
            # Pure solids and liquids use activity = 1: omit from Q, leave K unchanged
            # (e.g. Kw = [H+][OH-], Ksp = [Ca2+][F-]^2 with no H2O or CaF2 terms).
            A[:, j] = 0.0
        if getattr(compound, "excess", False):
            S[j, :] = 0.0
        for rxn in env.reactions:
            for species in rxn.reactants + rxn.products:
                if species["compound"] is compound and species.get("excess", False):
                    S[j, :] = 0.0
                    break

    for j in getattr(env, "buffer_indices", []):
        S[j, :] = 0.0

    from ._half_reaction import coupled_eh_mode, redox_reaction_indices

    redox_indices = np.array(redox_reaction_indices(env), dtype=int)
    couple_eh = coupled_eh_mode(env)

    ctx = EquilibriumContext(
        N=N, S=S, A=A, c0=c0, lnK=lnK, R=R, C=C,
        min_concentration=min_concentration,
        infinite_k_mask=infinite_k_mask,
        excess_source_mask=excess_source_mask,
        env=env,
        ln_gamma=np.zeros(C, dtype=float),
        redox_indices=redox_indices,
        couple_eh=couple_eh,
    )
    _update_activity(ctx, c0)
    return ctx


def _compute_lnQ(ctx: EquilibriumContext, c_safe: np.ndarray, c_actual: np.ndarray = None) -> np.ndarray:
    if c_actual is not None:
        _update_activity(ctx, c_actual)
    elif ctx.ln_gamma is None:
        _update_activity(ctx, c_safe)
    return ctx.A @ (np.log(c_safe) + ctx.ln_gamma)


def _compute_jacobian(ctx: EquilibriumContext, c_safe: np.ndarray, jacobian_scale: np.ndarray) -> np.ndarray:
    inv_c = 1.0 / c_safe
    J = ctx.A @ (inv_c[:, None] * ctx.S)
    return jacobian_scale[:, None] * J


def _reaction_quotient_errors(
    ctx: EquilibriumContext,
    x: np.ndarray,
) -> np.ndarray:
    """Per-reaction |Q/K - 1| for finite-K reactions (always computed, not masked)."""
    c = ctx.c0 + ctx.S @ x
    c_safe = _safe_concentrations(c, ctx.min_concentration)
    lnQ = _compute_lnQ(ctx, c_safe, c)
    q_over_k = np.exp(lnQ - ctx.lnK)
    errors = np.abs(q_over_k - 1.0)
    errors[ctx.infinite_k_mask] = 0.0
    return errors


def _finite_reaction_mask(ctx: EquilibriumContext) -> np.ndarray:
    return ~ctx.infinite_k_mask


def _reaction_is_negligible(
    ctx: EquilibriumContext,
    reaction_index: int,
    c_safe: np.ndarray,
    *,
    negligible_threshold: float = 1e-8,
) -> bool:
    if ctx.excess_source_mask[reaction_index]:
        return False
    active = np.abs(ctx.A[reaction_index, :]) > 0
    return active.any() and float(np.max(c_safe[active])) < negligible_threshold


def _active_residual(
    ctx: EquilibriumContext,
    c_safe: np.ndarray,
    loss_fn: EquilibriumLoss,
    *,
    c_actual: np.ndarray = None,
    negligible_threshold: float = 1e-8,
) -> np.ndarray:
    if c_actual is None:
        c_actual = c_safe
    lnQ = _compute_lnQ(ctx, c_safe, c_actual)
    residual = loss_fn.residual(lnQ, ctx.lnK)
    residual = residual.copy()
    for i in range(ctx.R):
        if ctx.infinite_k_mask[i] or _reaction_is_negligible(
            ctx, i, c_safe, negligible_threshold=negligible_threshold
        ):
            residual[i] = 0.0
    return residual


def _active_jacobian(
    ctx: EquilibriumContext,
    c_safe: np.ndarray,
    loss_fn: EquilibriumLoss,
    residual: np.ndarray,
) -> np.ndarray:
    scale = loss_fn.jacobian_scale(residual)
    J = _compute_jacobian(ctx, c_safe, scale)
    J = J.copy()
    for i in range(ctx.R):
        if ctx.infinite_k_mask[i]:
            J[i, :] = 0.0
    return J


def _use_residual_tolerance(quotient_error_limit: Optional[float]) -> bool:
    return quotient_error_limit is None


def _quotient_error_converged(
    errors: np.ndarray,
    quotient_error_limit: Optional[float],
    infinite_k_mask: np.ndarray,
) -> bool:
    if quotient_error_limit is None:
        return False
    finite_errors = errors[~infinite_k_mask] if infinite_k_mask.any() else errors
    if finite_errors.size == 0:
        return True
    return float(np.max(finite_errors)) <= quotient_error_limit


def _finalize(ctx: EquilibriumContext, x: np.ndarray):
    c_final = ctx.c0 + ctx.S @ x
    c_final = np.maximum(c_final, 0.0)
    return c_final.tolist(), x


def _needs_warm_start(
    ctx: EquilibriumContext,
    x: np.ndarray,
    loss_fn: EquilibriumLoss,
    tol: float,
) -> bool:
    """True when zero-product concentrations block first-order methods."""
    c = ctx.c0 + ctx.S @ x
    c_safe = _safe_concentrations(c, ctx.min_concentration)
    residual = _active_residual(ctx, c_safe, loss_fn, c_actual=c)
    if float(np.linalg.norm(residual, ord=2)) < max(10.0 * tol, 1e-6):
        return False
    zero_in_quotient = any(c[j] <= 0 and np.any(ctx.A[:, j] != 0) for j in range(ctx.C))
    return zero_in_quotient


def _newton_warm_start(
    ctx: EquilibriumContext,
    x: np.ndarray,
    loss_fn: EquilibriumLoss,
    *,
    max_steps: int = 8,
    learning_rate: float = 1.0,
    backtrack_beta: float = 0.5,
) -> np.ndarray:
    """A few damped Newton steps to escape zero-product singularities."""
    for _ in range(max_steps):
        c = ctx.c0 + ctx.S @ x
        c_safe = _safe_concentrations(c, ctx.min_concentration)
        residual = _active_residual(ctx, c_safe, loss_fn, c_actual=c)
        residual_norm = float(np.linalg.norm(residual, ord=2))
        if residual_norm < 1.0:
            break

        J = _active_jacobian(ctx, c_safe, loss_fn, residual)
        try:
            dx, *_ = np.linalg.lstsq(J, residual, rcond=None)
        except Exception:
            dx = np.linalg.pinv(J) @ residual
        dx[ctx.infinite_k_mask] = 0.0

        step = learning_rate
        improved = False
        while step >= 1e-12:
            x_new = x - step * dx
            x_new = _apply_infinite_k_extents(ctx, x_new)
            c_new = ctx.c0 + ctx.S @ x_new
            if np.all(c_new >= -1e-15):
                c_new_safe = _safe_concentrations(c_new, ctx.min_concentration)
                r_new = _active_residual(ctx, c_new_safe, loss_fn, c_actual=c_new)
                if float(np.linalg.norm(r_new, ord=2)) < residual_norm or step < 1e-10:
                    x = x_new
                    improved = True
                    break
            step *= backtrack_beta
        if not improved:
            break
    return x


def _prepare_first_order_extents(
    ctx: EquilibriumContext,
    x: np.ndarray,
    loss_fn: EquilibriumLoss,
    tol: float,
) -> np.ndarray:
    if _needs_warm_start(ctx, x, loss_fn, tol):
        x = _newton_warm_start(ctx, x, loss_fn)
    return x


def _build_result(
    env,
    ctx: EquilibriumContext,
    x: np.ndarray,
    concentrations: list[float],
    loss_fn: EquilibriumLoss,
    tol: float,
    stop_reason: str,
    iterations: int,
    quotient_error_limit: Optional[float],
) -> EquilibriumResult:
    errors = _reaction_quotient_errors(ctx, x)
    finite_errors = errors[_finite_reaction_mask(ctx)]
    max_error = float(np.max(finite_errors)) if finite_errors.size else 0.0

    c = ctx.c0 + ctx.S @ x
    c_safe = _safe_concentrations(c, ctx.min_concentration)
    residual = _active_residual(ctx, c_safe, loss_fn, c_actual=c)
    residual_norm = float(np.linalg.norm(residual, ord=2))
    lnQ = _compute_lnQ(ctx, c_safe, c)
    q_over_k = np.exp(lnQ - ctx.lnK)
    q_over_k = q_over_k.copy()
    q_over_k[ctx.infinite_k_mask] = 0.0

    if quotient_error_limit is not None:
        criterion_type = "quotient_error"
        criterion_value = max_error
        criterion_limit = float(quotient_error_limit)
        criterion_met = max_error <= quotient_error_limit
    else:
        criterion_type = "residual_tol"
        criterion_value = residual_norm
        criterion_limit = float(tol)
        criterion_met = residual_norm < tol

    if stop_reason == "quotient_error_limit":
        converged = True
    elif stop_reason == "residual_tol":
        converged = True
    else:
        converged = False

    return EquilibriumResult(
        concentrations=concentrations,
        compounds=[compound.formula for compound in env.compounds],
        reaction_extents=x.tolist(),
        reaction_quotient_error=errors.tolist(),
        max_reaction_quotient_error=max_error,
        reaction_quotient_ratio=q_over_k.tolist(),
        converged=converged,
        stop_reason=stop_reason,
        iterations=iterations,
        criterion_met=criterion_met,
        criterion_type=criterion_type,
        criterion_value=criterion_value,
        criterion_limit=criterion_limit,
        electrode_Eh=getattr(ctx, "solved_electrode_Eh", None),
    )


def _run_bgd(
    ctx: EquilibriumContext,
    loss_fn: EquilibriumLoss,
    max_iter: int,
    learning_rate: float,
    tol: float,
    backtrack_beta: float,
    quotient_error_limit: Optional[float],
):
    x = _initialize_extents(ctx)
    x = _prepare_first_order_extents(ctx, x, loss_fn, tol)
    stop_reason = "max_iter"

    for iteration in range(max_iter):
        c = ctx.c0 + ctx.S @ x
        c_safe = _safe_concentrations(c, ctx.min_concentration)

        lnQ = _compute_lnQ(ctx, c_safe, c)
        residual = _active_residual(ctx, c_safe, loss_fn, c_actual=c)

        errors = _reaction_quotient_errors(ctx, x)
        if _quotient_error_converged(errors, quotient_error_limit, ctx.infinite_k_mask):
            stop_reason = "quotient_error_limit"
            break

        if _use_residual_tolerance(quotient_error_limit) and np.linalg.norm(residual, ord=2) < tol:
            stop_reason = "residual_tol"
            break

        J = _active_jacobian(ctx, c_safe, loss_fn, residual)
        grad = J.T @ loss_fn.grad_weights(residual)
        grad[ctx.infinite_k_mask] = 0.0

        if _use_residual_tolerance(quotient_error_limit) and np.linalg.norm(grad, ord=2) < tol:
            stop_reason = "residual_tol"
            break

        step = learning_rate
        f_curr = loss_fn.objective(residual)
        while True:
            x_new = x - step * grad
            x_new = _apply_infinite_k_extents(ctx, x_new)
            c_new = ctx.c0 + ctx.S @ x_new
            if np.all(c_new >= -1e-15):
                c_new_safe = _safe_concentrations(c_new, ctx.min_concentration)
                r_new = _active_residual(ctx, c_new_safe, loss_fn, c_actual=c_new)
                f_new = loss_fn.objective(r_new)
                if f_new <= f_curr or step < 1e-12:
                    x = x_new
                    break
            step *= backtrack_beta
    else:
        iteration = max_iter - 1

    concentrations, x = _finalize(ctx, x)
    return concentrations, x, stop_reason, iteration + 1


def _run_sgd(
    ctx: EquilibriumContext,
    loss_fn: EquilibriumLoss,
    max_iter: int,
    learning_rate: float,
    tol: float,
    backtrack_beta: float,
    quotient_error_limit: Optional[float],
):
    x = _initialize_extents(ctx)
    x = _prepare_first_order_extents(ctx, x, loss_fn, tol)
    stop_reason = "max_iter"

    for iteration in range(max_iter):
        order = np.random.permutation(ctx.R)
        any_update = False

        for i in order:
            if ctx.infinite_k_mask[i]:
                continue
            c = ctx.c0 + ctx.S @ x
            c_safe = _safe_concentrations(c, ctx.min_concentration)
            inv_c = 1.0 / c_safe

            _update_activity(ctx, c)
            lnQ_i = ctx.A[i, :] @ (np.log(c_safe) + ctx.ln_gamma)
            r_i = loss_fn.residual(np.array([lnQ_i]), np.array([ctx.lnK[i]]))[0]
            if _use_residual_tolerance(quotient_error_limit) and abs(r_i) < tol:
                continue

            scale_i = loss_fn.jacobian_scale(np.array([r_i]))[0]
            J_i = scale_i * (ctx.A[i, :] @ (inv_c[:, None] * ctx.S))
            grad_i = J_i * loss_fn.grad_weights(np.array([r_i]))[0]

            step = learning_rate
            f_curr = loss_fn.objective(np.array([r_i]))

            while True:
                x_new = x - step * grad_i
                x_new = _apply_infinite_k_extents(ctx, x_new)
                c_new = ctx.c0 + ctx.S @ x_new
                if np.all(c_new >= -1e-15):
                    c_new_safe = _safe_concentrations(c_new, ctx.min_concentration)
                    _update_activity(ctx, c_new)
                    lnQ_i_new = ctx.A[i, :] @ (np.log(c_new_safe) + ctx.ln_gamma)
                    r_i_new = loss_fn.residual(np.array([lnQ_i_new]), np.array([ctx.lnK[i]]))[0]
                    f_new = loss_fn.objective(np.array([r_i_new]))
                    if f_new <= f_curr or step < 1e-12:
                        x = x_new
                        any_update = True
                        break
                step *= backtrack_beta

        errors = _reaction_quotient_errors(ctx, x)
        if _quotient_error_converged(errors, quotient_error_limit, ctx.infinite_k_mask):
            stop_reason = "quotient_error_limit"
            break

        c_full = ctx.c0 + ctx.S @ x
        c_full_safe = _safe_concentrations(c_full, ctx.min_concentration)
        full_residual = _active_residual(ctx, c_full_safe, loss_fn, c_actual=c_full)

        if _use_residual_tolerance(quotient_error_limit) and np.linalg.norm(full_residual, ord=2) < tol:
            stop_reason = "residual_tol"
            break
        if not any_update:
            break
    else:
        iteration = max_iter - 1

    concentrations, x = _finalize(ctx, x)
    return concentrations, x, stop_reason, iteration + 1


def _run_newton(
    ctx: EquilibriumContext,
    loss_fn: EquilibriumLoss,
    max_iter: int,
    learning_rate: float,
    tol: float,
    backtrack_beta: float,
    quotient_error_limit: Optional[float],
):
    x = _initialize_extents(ctx)
    stop_reason = "max_iter"

    for iteration in range(max_iter):
        c = ctx.c0 + ctx.S @ x
        c_safe = _safe_concentrations(c, ctx.min_concentration)

        residual = _active_residual(ctx, c_safe, loss_fn, c_actual=c)

        errors = _reaction_quotient_errors(ctx, x)
        if _quotient_error_converged(errors, quotient_error_limit, ctx.infinite_k_mask):
            stop_reason = "quotient_error_limit"
            break

        if _use_residual_tolerance(quotient_error_limit) and np.linalg.norm(residual, ord=2) < tol:
            stop_reason = "residual_tol"
            break

        J = _active_jacobian(ctx, c_safe, loss_fn, residual)

        try:
            dx, *_ = np.linalg.lstsq(J, residual, rcond=None)
        except Exception:
            dx = np.linalg.pinv(J) @ residual
        dx[ctx.infinite_k_mask] = 0.0

        step = learning_rate
        f_curr = loss_fn.objective(residual)
        while True:
            x_new = x - step * dx
            x_new = _apply_infinite_k_extents(ctx, x_new)
            c_new = ctx.c0 + ctx.S @ x_new
            if np.all(c_new >= -1e-15):
                c_new_safe = _safe_concentrations(c_new, ctx.min_concentration)
                r_new = _active_residual(ctx, c_new_safe, loss_fn, c_actual=c_new)
                f_new = loss_fn.objective(r_new)
                if f_new <= f_curr or step < 1e-12:
                    x = x_new
                    break
            step *= backtrack_beta
    else:
        iteration = max_iter - 1

    concentrations, x = _finalize(ctx, x)
    return concentrations, x, stop_reason, iteration + 1


def _refresh_redox_lnK(ctx: EquilibriumContext, x: np.ndarray, Eh: float) -> None:
    from ._half_reaction import compute_pH, redox_lnK

    env = ctx.env
    c = ctx.c0 + ctx.S @ x
    pH = compute_pH(env, c)
    for hr in getattr(env, "half_reactions", None) or []:
        idx = hr._reaction_index
        if idx is None:
            continue
        ctx.lnK[idx] = redox_lnK(hr, Eh, pH, env.T)


def _coupled_eh_residual(
    ctx: EquilibriumContext,
    x: np.ndarray,
    Eh: float,
    loss_fn: EquilibriumLoss,
) -> np.ndarray:
    _refresh_redox_lnK(ctx, x, Eh)
    c = ctx.c0 + ctx.S @ x
    c_safe = _safe_concentrations(c, ctx.min_concentration)
    return _active_residual(ctx, c_safe, loss_fn, c_actual=c)


def _coupled_eh_jacobian(
    ctx: EquilibriumContext,
    x: np.ndarray,
    Eh: float,
    loss_fn: EquilibriumLoss,
    residual: np.ndarray,
) -> np.ndarray:
    from ._half_reaction import d_redox_lnK_dEh

    c = ctx.c0 + ctx.S @ x
    c_safe = _safe_concentrations(c, ctx.min_concentration)
    jx = _active_jacobian(ctx, c_safe, loss_fn, residual)
    deh = np.zeros(ctx.R, dtype=float)
    env = ctx.env
    for hr in getattr(env, "half_reactions", None) or []:
        idx = hr._reaction_index
        if idx is None or ctx.infinite_k_mask[idx]:
            continue
        deh[idx] = -d_redox_lnK_dEh(hr, env.T) * loss_fn.jacobian_scale(residual)[idx]
    return np.column_stack([jx, deh])


def _run_newton_coupled_eh(
    ctx: EquilibriumContext,
    loss_fn: EquilibriumLoss,
    max_iter: int,
    learning_rate: float,
    tol: float,
    backtrack_beta: float,
    quotient_error_limit: Optional[float],
):
    from ._half_reaction import initial_eh_guess

    x = _initialize_extents(ctx)
    Eh = initial_eh_guess(ctx.env)
    stop_reason = "max_iter"

    for iteration in range(max_iter):
        residual = _coupled_eh_residual(ctx, x, Eh, loss_fn)
        errors = _reaction_quotient_errors(ctx, x)
        if _quotient_error_converged(errors, quotient_error_limit, ctx.infinite_k_mask):
            stop_reason = "quotient_error_limit"
            break
        if _use_residual_tolerance(quotient_error_limit) and np.linalg.norm(residual, ord=2) < tol:
            stop_reason = "residual_tol"
            break

        J = _coupled_eh_jacobian(ctx, x, Eh, loss_fn, residual)
        try:
            step_vec, *_ = np.linalg.lstsq(J, residual, rcond=None)
        except Exception:
            step_vec = np.linalg.pinv(J) @ residual
        dx = step_vec[: ctx.R]
        dEh = float(step_vec[ctx.R]) if len(step_vec) > ctx.R else 0.0
        dx[ctx.infinite_k_mask] = 0.0

        step = learning_rate
        f_curr = loss_fn.objective(residual)
        while True:
            x_new = x - step * dx
            x_new = _apply_infinite_k_extents(ctx, x_new)
            Eh_new = Eh - step * dEh
            c_new = ctx.c0 + ctx.S @ x_new
            if np.all(c_new >= -1e-15):
                c_new_safe = _safe_concentrations(c_new, ctx.min_concentration)
                r_new = _coupled_eh_residual(ctx, x_new, Eh_new, loss_fn)
                f_new = loss_fn.objective(r_new)
                if f_new <= f_curr or step < 1e-12:
                    x = x_new
                    Eh = Eh_new
                    break
            step *= backtrack_beta
    else:
        iteration = max_iter - 1

    ctx.solved_electrode_Eh = float(Eh)
    concentrations, x = _finalize(ctx, x)
    return concentrations, x, stop_reason, iteration + 1


def solve_equilibrium(
    env,
    *,
    method: str = "bgd",
    loss: str = "log_quotient",
    max_iter: Optional[int] = None,
    learning_rate: Optional[float] = None,
    tol: Optional[float] = None,
    backtrack_beta: float = 0.5,
    min_concentration: float = 1e-12,
    quotient_error_limit: Optional[float] = None,
    huber_delta: float = 1.0,
    return_details: bool = False,
) -> Union[list[float], EquilibriumResult]:
    if method not in VALID_METHODS:
        raise ValueError(f"Invalid method {method!r}. Choose from {VALID_METHODS}.")
    if loss not in LOSS_REGISTRY:
        raise ValueError(f"Invalid loss {loss!r}. Choose from {tuple(LOSS_REGISTRY.keys())}.")

    defaults = METHOD_DEFAULTS[method]
    max_iter = defaults["max_iter"] if max_iter is None else max_iter
    learning_rate = defaults["learning_rate"] if learning_rate is None else learning_rate
    tol = defaults["tol"] if tol is None else tol

    from ._half_reaction import apply_electrode_potential, coupled_eh_mode

    if not coupled_eh_mode(env):
        apply_electrode_potential(env)

    ctx = _build_context(env, min_concentration)
    loss_fn = EquilibriumLoss(loss, delta=huber_delta)

    solver_kwargs = {
        "max_iter": max_iter,
        "learning_rate": learning_rate,
        "tol": tol,
        "backtrack_beta": backtrack_beta,
        "quotient_error_limit": quotient_error_limit,
    }

    solved_eh = getattr(env, "electrode_Eh", None)
    if ctx.couple_eh:
        concentrations, x, stop_reason, iterations = _run_newton_coupled_eh(ctx, loss_fn, **solver_kwargs)
        solved_eh = getattr(ctx, "solved_electrode_Eh", solved_eh)
    elif method == "bgd":
        concentrations, x, stop_reason, iterations = _run_bgd(ctx, loss_fn, **solver_kwargs)
    elif method == "sgd":
        concentrations, x, stop_reason, iterations = _run_sgd(ctx, loss_fn, **solver_kwargs)
    else:
        concentrations, x, stop_reason, iterations = _run_newton(ctx, loss_fn, **solver_kwargs)

    if solved_eh is None and getattr(env, "electrode_Eh", None) is not None:
        solved_eh = env.electrode_Eh

    env._equilibrium_x_solution = x

    result = _build_result(
        env,
        ctx,
        x,
        concentrations,
        loss_fn,
        tol,
        stop_reason,
        iterations,
        quotient_error_limit,
    )
    result.electrode_Eh = solved_eh

    if return_details:
        return result
    return concentrations
