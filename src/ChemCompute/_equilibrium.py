from dataclasses import dataclass
from typing import Optional

import numpy as np

METHOD_DEFAULTS = {
    "bgd": {"max_iter": 5000, "learning_rate": 0.1, "tol": 1e-8},
    "sgd": {"max_iter": 5000, "learning_rate": 0.1, "tol": 1e-8},
    "newton": {"max_iter": 200, "learning_rate": 1.0, "tol": 1e-10},
}

VALID_METHODS = tuple(METHOD_DEFAULTS.keys())


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


def _concentration_converged(
    c_new: np.ndarray,
    c_prev: Optional[np.ndarray],
    min_concentration: float,
    concentration_error_limit: Optional[float],
) -> bool:
    if c_prev is None or concentration_error_limit is None:
        return False
    denom = np.maximum(np.abs(c_prev), min_concentration)
    rel_change = np.max(np.abs(c_new - c_prev) / denom)
    return rel_change < concentration_error_limit


def _finalize(ctx: EquilibriumContext, x: np.ndarray):
    c_final = ctx.c0 + ctx.S @ x
    c_final = np.maximum(c_final, 0.0)
    return c_final.tolist(), x


def _run_bgd(
    ctx: EquilibriumContext,
    loss_fn: EquilibriumLoss,
    max_iter: int,
    learning_rate: float,
    tol: float,
    backtrack_beta: float,
    concentration_error_limit: Optional[float],
):
    x = np.zeros(ctx.R, dtype=float)
    c_prev = None

    for _ in range(max_iter):
        c = ctx.c0 + ctx.S @ x
        c_safe = np.maximum(c, ctx.min_concentration)

        lnQ = _compute_lnQ(ctx, c_safe)
        residual = loss_fn.residual(lnQ, ctx.lnK)

        if np.linalg.norm(residual, ord=2) < tol:
            break
        if _concentration_converged(c, c_prev, ctx.min_concentration, concentration_error_limit):
            break

        scale = loss_fn.jacobian_scale(residual)
        J = _compute_jacobian(ctx, c_safe, scale)
        grad = J.T @ loss_fn.grad_weights(residual)

        if np.linalg.norm(grad, ord=2) < tol:
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
                    c_prev = c.copy()
                    x = x_new
                    break
            step *= backtrack_beta

    return _finalize(ctx, x)


def _run_sgd(
    ctx: EquilibriumContext,
    loss_fn: EquilibriumLoss,
    max_iter: int,
    learning_rate: float,
    tol: float,
    backtrack_beta: float,
    concentration_error_limit: Optional[float],
):
    x = np.zeros(ctx.R, dtype=float)
    c_prev = None

    for _ in range(max_iter):
        order = np.random.permutation(ctx.R)
        any_update = False

        for i in order:
            c = ctx.c0 + ctx.S @ x
            c_safe = np.maximum(c, ctx.min_concentration)
            inv_c = 1.0 / c_safe

            lnQ_i = ctx.A[i, :] @ np.log(c_safe)
            r_i = loss_fn.residual(np.array([lnQ_i]), np.array([ctx.lnK[i]]))[0]
            if abs(r_i) < tol:
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
                        c_prev = c.copy()
                        x = x_new
                        any_update = True
                        break
                step *= backtrack_beta

        c_full = ctx.c0 + ctx.S @ x
        c_full_safe = np.maximum(c_full, ctx.min_concentration)
        full_lnQ = _compute_lnQ(ctx, c_full_safe)
        full_residual = loss_fn.residual(full_lnQ, ctx.lnK)

        if np.linalg.norm(full_residual, ord=2) < tol:
            break
        if _concentration_converged(c_full, c_prev, ctx.min_concentration, concentration_error_limit):
            break
        if not any_update:
            break

    return _finalize(ctx, x)


def _run_newton(
    ctx: EquilibriumContext,
    loss_fn: EquilibriumLoss,
    max_iter: int,
    learning_rate: float,
    tol: float,
    backtrack_beta: float,
    concentration_error_limit: Optional[float],
):
    x = np.zeros(ctx.R, dtype=float)
    c_prev = None

    for _ in range(max_iter):
        c = ctx.c0 + ctx.S @ x
        c_safe = np.maximum(c, ctx.min_concentration)

        lnQ = _compute_lnQ(ctx, c_safe)
        residual = loss_fn.residual(lnQ, ctx.lnK)

        if np.linalg.norm(residual, ord=2) < tol:
            break
        if _concentration_converged(c, c_prev, ctx.min_concentration, concentration_error_limit):
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
                    c_prev = c.copy()
                    x = x_new
                    break
            step *= backtrack_beta

    return _finalize(ctx, x)


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
    concentration_error_limit: Optional[float] = None,
    huber_delta: float = 1.0,
) -> list[float]:
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
        "concentration_error_limit": concentration_error_limit,
    }

    if method == "bgd":
        concentrations, x = _run_bgd(ctx, loss_fn, **solver_kwargs)
    elif method == "sgd":
        concentrations, x = _run_sgd(ctx, loss_fn, **solver_kwargs)
    else:
        concentrations, x = _run_newton(ctx, loss_fn, **solver_kwargs)

    env._equilibrium_x_solution = x
    return concentrations
