from dataclasses import dataclass
from typing import Optional, Union

import numpy as np

METHOD_DEFAULTS = {
    "bgd": {"max_iter": 5000, "learning_rate": 0.1, "tol": 1e-8},
    "sgd": {"max_iter": 5000, "learning_rate": 0.1, "tol": 1e-8},
    "newton": {"max_iter": 200, "learning_rate": 1.0, "tol": 1e-10},
}

VALID_METHODS = tuple(METHOD_DEFAULTS.keys())
EXTENT_FLOOR = 1e-12
JACOBIAN_DIAG_FLOOR = 1e-30


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
    concentrations: list[float]
    compounds: list[str]
    reaction_extents: list[float]
    reaction_extent_percent: list[float]
    max_reaction_extent_percent: float
    reaction_quotient_ratio: list[float]
    converged: bool
    stop_reason: str
    iterations: int

    @property
    def q_over_k(self) -> list[float]:
        """Per-reaction Q/K ratio at the solution (1.0 means equilibrium)."""
        return self.reaction_quotient_ratio


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


def _build_context(env, min_concentration: float) -> EquilibriumContext:
    N = env.stoichiometric_coefficient_array
    S = N.T
    R, C = N.shape

    c0 = np.array(env.concentrations, dtype=float)
    A = -N.astype(float)

    for j, compound in enumerate(env.compounds):
        ph = compound.phase(env.T)
        if ph in ("s", "l"):
            A[:, j] = 0.0

    K_vec = np.array([max(rxn.K, 1e-300) for rxn in env.reactions], dtype=float)
    lnK = np.log(K_vec)

    return EquilibriumContext(N=N, S=S, A=A, c0=c0, lnK=lnK, R=R, C=C, min_concentration=min_concentration)


def _compute_lnQ(ctx: EquilibriumContext, c_safe: np.ndarray) -> np.ndarray:
    return ctx.A @ np.log(c_safe)


def _compute_jacobian(ctx: EquilibriumContext, c_safe: np.ndarray, jacobian_scale: np.ndarray) -> np.ndarray:
    inv_c = 1.0 / c_safe
    J = ctx.A @ (inv_c[:, None] * ctx.S)
    return jacobian_scale[:, None] * J


def _residual_i_at_delta(
    ctx: EquilibriumContext,
    x: np.ndarray,
    reaction_index: int,
    delta: float,
    loss_fn: EquilibriumLoss,
) -> float:
    x_try = x.copy()
    x_try[reaction_index] += delta
    c = ctx.c0 + ctx.S @ x_try
    c_safe = np.maximum(c, ctx.min_concentration)
    lnQ = _compute_lnQ(ctx, c_safe)
    residual = loss_fn.residual(lnQ, ctx.lnK)
    return float(residual[reaction_index])


def _bisection_delta(
    ctx: EquilibriumContext,
    x: np.ndarray,
    reaction_index: int,
    loss_fn: EquilibriumLoss,
    tol: float,
) -> float:
    f0 = _residual_i_at_delta(ctx, x, reaction_index, 0.0, loss_fn)
    if abs(f0) < tol:
        return 0.0

    low, high = -1.0, 1.0
    f_low = _residual_i_at_delta(ctx, x, reaction_index, low, loss_fn)
    f_high = _residual_i_at_delta(ctx, x, reaction_index, high, loss_fn)

    expand = 1.0
    while f_low * f_high > 0 and expand < 1e6:
        expand *= 2.0
        low = -expand
        high = expand
        f_low = _residual_i_at_delta(ctx, x, reaction_index, low, loss_fn)
        f_high = _residual_i_at_delta(ctx, x, reaction_index, high, loss_fn)

    if f_low * f_high > 0:
        return 0.0

    for _ in range(80):
        mid = 0.5 * (low + high)
        f_mid = _residual_i_at_delta(ctx, x, reaction_index, mid, loss_fn)
        if abs(f_mid) < tol or (high - low) < 1e-14:
            return mid
        if f_low * f_mid <= 0:
            high = mid
            f_high = f_mid
        else:
            low = mid
            f_low = f_mid
    return 0.5 * (low + high)


def _reaction_extent_gaps(
    ctx: EquilibriumContext,
    x: np.ndarray,
    loss_fn: EquilibriumLoss,
    tol: float,
) -> tuple[np.ndarray, np.ndarray]:
    c = ctx.c0 + ctx.S @ x
    c_safe = np.maximum(c, ctx.min_concentration)
    lnQ = _compute_lnQ(ctx, c_safe)
    residual = loss_fn.residual(lnQ, ctx.lnK)
    scale = loss_fn.jacobian_scale(residual)
    J = _compute_jacobian(ctx, c_safe, scale)

    delta_x = np.zeros(ctx.R, dtype=float)
    for i in range(ctx.R):
        r_i = residual[i]
        if abs(r_i) < tol:
            continue

        j_ii = J[i, i]
        if abs(j_ii) > JACOBIAN_DIAG_FLOOR:
            delta_x[i] = -r_i / j_ii
        else:
            delta_x[i] = _bisection_delta(ctx, x, i, loss_fn, tol)

    denom = np.maximum(np.abs(x), EXTENT_FLOOR)
    percent = np.abs(delta_x) / denom
    return delta_x, percent


def _use_residual_tolerance(reaction_extent_error_limit: Optional[float]) -> bool:
    return reaction_extent_error_limit is None


def _reaction_extent_converged(
    percent: np.ndarray,
    reaction_extent_error_limit: Optional[float],
) -> bool:
    if reaction_extent_error_limit is None:
        return False
    return float(np.max(percent)) <= reaction_extent_error_limit


def _finalize(ctx: EquilibriumContext, x: np.ndarray):
    c_final = ctx.c0 + ctx.S @ x
    c_final = np.maximum(c_final, 0.0)
    return c_final.tolist(), x


def _build_result(
    env,
    ctx: EquilibriumContext,
    x: np.ndarray,
    concentrations: list[float],
    loss_fn: EquilibriumLoss,
    tol: float,
    stop_reason: str,
    iterations: int,
    reaction_extent_error_limit: Optional[float],
) -> EquilibriumResult:
    _, percent = _reaction_extent_gaps(ctx, x, loss_fn, tol)
    max_percent = float(np.max(percent)) if percent.size else 0.0

    c = ctx.c0 + ctx.S @ x
    c_safe = np.maximum(c, ctx.min_concentration)
    lnQ = _compute_lnQ(ctx, c_safe)
    q_over_k = np.exp(lnQ - ctx.lnK)

    if stop_reason == "reaction_extent_limit":
        converged = True
    elif stop_reason == "residual_tol":
        converged = True
    else:
        converged = False

    return EquilibriumResult(
        concentrations=concentrations,
        compounds=[compound.formula for compound in env.compounds],
        reaction_extents=x.tolist(),
        reaction_extent_percent=percent.tolist(),
        max_reaction_extent_percent=max_percent,
        reaction_quotient_ratio=q_over_k.tolist(),
        converged=converged,
        stop_reason=stop_reason,
        iterations=iterations,
    )


def _run_bgd(
    ctx: EquilibriumContext,
    loss_fn: EquilibriumLoss,
    max_iter: int,
    learning_rate: float,
    tol: float,
    backtrack_beta: float,
    reaction_extent_error_limit: Optional[float],
):
    x = np.zeros(ctx.R, dtype=float)
    stop_reason = "max_iter"

    for iteration in range(max_iter):
        c = ctx.c0 + ctx.S @ x
        c_safe = np.maximum(c, ctx.min_concentration)

        lnQ = _compute_lnQ(ctx, c_safe)
        residual = loss_fn.residual(lnQ, ctx.lnK)

        _, percent = _reaction_extent_gaps(ctx, x, loss_fn, tol)
        if _reaction_extent_converged(percent, reaction_extent_error_limit):
            stop_reason = "reaction_extent_limit"
            break

        if _use_residual_tolerance(reaction_extent_error_limit) and np.linalg.norm(residual, ord=2) < tol:
            stop_reason = "residual_tol"
            break

        scale = loss_fn.jacobian_scale(residual)
        J = _compute_jacobian(ctx, c_safe, scale)
        grad = J.T @ loss_fn.grad_weights(residual)

        if _use_residual_tolerance(reaction_extent_error_limit) and np.linalg.norm(grad, ord=2) < tol:
            stop_reason = "residual_tol"
            break

        step = learning_rate
        f_curr = loss_fn.objective(residual)
        while True:
            x_new = x - step * grad
            c_new = ctx.c0 + ctx.S @ x_new
            if np.all(c_new >= -1e-15):
                c_new_safe = np.maximum(c_new, ctx.min_concentration)
                lnQ_new = _compute_lnQ(ctx, c_new_safe)
                r_new = loss_fn.residual(lnQ_new, ctx.lnK)
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
    reaction_extent_error_limit: Optional[float],
):
    x = np.zeros(ctx.R, dtype=float)
    stop_reason = "max_iter"

    for iteration in range(max_iter):
        order = np.random.permutation(ctx.R)
        any_update = False

        for i in order:
            c = ctx.c0 + ctx.S @ x
            c_safe = np.maximum(c, ctx.min_concentration)
            inv_c = 1.0 / c_safe

            lnQ_i = ctx.A[i, :] @ np.log(c_safe)
            r_i = loss_fn.residual(np.array([lnQ_i]), np.array([ctx.lnK[i]]))[0]
            if _use_residual_tolerance(reaction_extent_error_limit) and abs(r_i) < tol:
                continue

            scale_i = loss_fn.jacobian_scale(np.array([r_i]))[0]
            J_i = scale_i * (ctx.A[i, :] @ (inv_c[:, None] * ctx.S))
            grad_i = J_i * loss_fn.grad_weights(np.array([r_i]))[0]

            step = learning_rate
            f_curr = loss_fn.objective(np.array([r_i]))

            while True:
                x_new = x - step * grad_i
                c_new = ctx.c0 + ctx.S @ x_new
                if np.all(c_new >= -1e-15):
                    c_new_safe = np.maximum(c_new, ctx.min_concentration)
                    lnQ_i_new = ctx.A[i, :] @ np.log(c_new_safe)
                    r_i_new = loss_fn.residual(np.array([lnQ_i_new]), np.array([ctx.lnK[i]]))[0]
                    f_new = loss_fn.objective(np.array([r_i_new]))
                    if f_new <= f_curr or step < 1e-12:
                        x = x_new
                        any_update = True
                        break
                step *= backtrack_beta

        _, percent = _reaction_extent_gaps(ctx, x, loss_fn, tol)
        if _reaction_extent_converged(percent, reaction_extent_error_limit):
            stop_reason = "reaction_extent_limit"
            break

        c_full = ctx.c0 + ctx.S @ x
        c_full_safe = np.maximum(c_full, ctx.min_concentration)
        full_lnQ = _compute_lnQ(ctx, c_full_safe)
        full_residual = loss_fn.residual(full_lnQ, ctx.lnK)

        if _use_residual_tolerance(reaction_extent_error_limit) and np.linalg.norm(full_residual, ord=2) < tol:
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
    reaction_extent_error_limit: Optional[float],
):
    x = np.zeros(ctx.R, dtype=float)
    stop_reason = "max_iter"

    for iteration in range(max_iter):
        c = ctx.c0 + ctx.S @ x
        c_safe = np.maximum(c, ctx.min_concentration)

        lnQ = _compute_lnQ(ctx, c_safe)
        residual = loss_fn.residual(lnQ, ctx.lnK)

        _, percent = _reaction_extent_gaps(ctx, x, loss_fn, tol)
        if _reaction_extent_converged(percent, reaction_extent_error_limit):
            stop_reason = "reaction_extent_limit"
            break

        if _use_residual_tolerance(reaction_extent_error_limit) and np.linalg.norm(residual, ord=2) < tol:
            stop_reason = "residual_tol"
            break

        scale = loss_fn.jacobian_scale(residual)
        J = _compute_jacobian(ctx, c_safe, scale)

        try:
            dx, *_ = np.linalg.lstsq(J, residual, rcond=None)
        except Exception:
            dx = np.linalg.pinv(J) @ residual

        step = learning_rate
        f_curr = loss_fn.objective(residual)
        while True:
            x_new = x - step * dx
            c_new = ctx.c0 + ctx.S @ x_new
            if np.all(c_new >= -1e-15):
                c_new_safe = np.maximum(c_new, ctx.min_concentration)
                lnQ_new = _compute_lnQ(ctx, c_new_safe)
                r_new = loss_fn.residual(lnQ_new, ctx.lnK)
                f_new = loss_fn.objective(r_new)
                if f_new <= f_curr or step < 1e-12:
                    x = x_new
                    break
            step *= backtrack_beta
    else:
        iteration = max_iter - 1

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
    reaction_extent_error_limit: Optional[float] = None,
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

    ctx = _build_context(env, min_concentration)
    loss_fn = EquilibriumLoss(loss, delta=huber_delta)

    solver_kwargs = {
        "max_iter": max_iter,
        "learning_rate": learning_rate,
        "tol": tol,
        "backtrack_beta": backtrack_beta,
        "reaction_extent_error_limit": reaction_extent_error_limit,
    }

    if method == "bgd":
        concentrations, x, stop_reason, iterations = _run_bgd(ctx, loss_fn, **solver_kwargs)
    elif method == "sgd":
        concentrations, x, stop_reason, iterations = _run_sgd(ctx, loss_fn, **solver_kwargs)
    else:
        concentrations, x, stop_reason, iterations = _run_newton(ctx, loss_fn, **solver_kwargs)

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
        reaction_extent_error_limit,
    )

    if return_details:
        return result
    return concentrations
