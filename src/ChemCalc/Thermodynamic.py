from _general import Enviroment
import numpy as np


class EquilibriumCalculator:
    def __init__(self, method_of_calculation: str = "bgd"):
        self.method_of_calculation = method_of_calculation
        self.fitted = False
    def _generate_concentration_equations(self):
        # Start with a copy of the current concentrations as strings
        concentration_eq = [str(value) for value in self.env.concentrations]

        # For each reaction, adjust each compound's expression by its stoichiometric change
        for r_index, reaction in enumerate(self.env.reactions, start=1):
            for idx, compound in enumerate(self.env.compounds):
                coeff = 0
                # Check reactants (positive coefficient)
                for reactant in reaction.reactants:
                    if reactant["compound"] == compound:
                        coeff += reactant["stoichiometric_coefficient"]
                        break
                # Check products (negative coefficient)
                for product in reaction.products:
                    if product["compound"] == compound:
                        coeff -= product["stoichiometric_coefficient"]
                        break
                if coeff != 0:
                    concentration_eq[idx] += f" + ({coeff}x{r_index})"
        return concentration_eq

    def fit(self, env: Enviroment):
        self.env = env
        self.concentration_equation = self._generate_concentration_equations()
        self.fitted = True
    def calculate(self,
                  max_iter: int = 5000,
                learning_rate: float = 0.1,
                tol: float = 1e-8,
                backtrack_beta: float = 0.5,
                min_concentration: float = 1e-12):
        if self.fitted == False:
            raise ValueError("Environment not fitted")
        if self.method_of_calculation == "bgd":
            return self._caculate_by_batch_gradient_descent(max_iter, learning_rate, tol, backtrack_beta, min_concentration)
        elif self.method_of_calculation == "sgd":
            return self._caculate_by_stochastic_gradient_descent(max_iter, learning_rate, tol, backtrack_beta, min_concentration)
        elif self.method_of_calculation == "newton":
            return self._calculate_by_newton(max_iter, learning_rate, tol, backtrack_beta, min_concentration)
        return None

    def _caculate_by_batch_gradient_descent(self,
                                            max_iter: int = 5000,
                                            learning_rate: float = 0.1,
                                            tol: float = 1e-8,
                                            backtrack_beta: float = 0.5,
                                            min_concentration: float = 1e-12):
        
        N = self.env.stoichiometric_coefficient_array  
        S = N.T 
        R = N.shape[0]
        C = N.shape[1]

        # Initial concentrations vector
        c0 = np.array(self.env.concentrations, dtype=float) 

        # Build mass-action exponent matrix A: products +, reactants -
        A = -N.astype(float)

        # Phase handling: exclude s/l from equilibrium expression by zeroing their columns in A
        phase_include_mask = np.ones(C, dtype=bool)
        for j, compound in enumerate(self.env.compounds):
            ph = compound.phase(self.env.T)
            if ph in ("s", "l"):
                phase_include_mask[j] = False
        for j in range(C):
            if not phase_include_mask[j]:
                A[:, j] = 0.0

        # Equilibrium constants vector
        K_vec = np.array([max(rxn.K, 1e-300) for rxn in self.env.reactions], dtype=float)
        lnK = np.log(K_vec)

        # Optimize extents x
        x = np.zeros(R, dtype=float)

        for _ in range(max_iter):
            c = c0 + S @ x  # (C,)
            # Ensure strictly positive for log; use floor at min_concentration for stability
            c_safe = np.maximum(c, min_concentration)

            lnQ = A @ np.log(c_safe) 
            residual = lnQ - lnK

            # Check residual convergence
            if np.linalg.norm(residual, ord=2) < tol:
                break

            # Jacobian J = A @ diag(1/c) @ S  (R x R)
            inv_c = 1.0 / c_safe
            J = A @ (inv_c[:, None] * S)  # broadcasting builds A @ diag(inv_c) @ S
            grad = J.T @ residual  # (R,)

            # Gradient convergence
            if np.linalg.norm(grad, ord=2) < tol:
                break

            # Take step with backtracking to preserve non-negativity and reduce objective
            step = learning_rate
            f_curr = 0.5 * np.dot(residual, residual)
            while True:
                x_new = x - step * grad
                c_new = c0 + S @ x_new
                if np.all(c_new >= -1e-15):  # allow tiny numerical negative, will be floored for logs
                    c_new_safe = np.maximum(c_new, min_concentration)
                    lnQ_new = A @ np.log(c_new_safe)
                    r_new = lnQ_new - lnK
                    f_new = 0.5 * np.dot(r_new, r_new)
                    if f_new <= f_curr or step < 1e-12:
                        x = x_new
                        break
                step *= backtrack_beta

        # Final concentrations
        c_final = c0 + S @ x
        c_final = np.maximum(c_final, 0.0)

        # Save state and return result
        self.x_solution = x
        self.fitted = True
        # Return as list aligned with env.compounds
        return c_final.tolist()

    def _calculate_by_stochastic_gradient_descent(self,
                                                  max_iter: int = 5000,
                                                  learning_rate: float = 0.1,
                                                  tol: float = 1e-8,
                                                  backtrack_beta: float = 0.5,
                                                  min_concentration: float = 1e-12):
        # Stoichiometry and mappings
        N = self.env.stoichiometric_coefficient_array
        S = N.T
        R = N.shape[0]
        C = N.shape[1]

        c0 = np.array(self.env.concentrations, dtype=float)

        # Mass-action exponents, exclude s/l phases
        A = -N.astype(float)
        phase_include_mask = np.ones(C, dtype=bool)
        for j, compound in enumerate(self.env.compounds):
            ph = compound.phase(self.env.T)
            if ph in ("s", "l"):
                phase_include_mask[j] = False
        for j in range(C):
            if not phase_include_mask[j]:
                A[:, j] = 0.0

        K_vec = np.array([max(rxn.K, 1e-300) for rxn in self.env.reactions], dtype=float)
        lnK = np.log(K_vec)

        x = np.zeros(R, dtype=float)

        for _ in range(max_iter):
            order = np.random.permutation(R)
            any_update = False
            for i in order:
                c = c0 + S @ x
                c_safe = np.maximum(c, min_concentration)
                inv_c = 1.0 / c_safe

                a_i = A[i, :]
                lnQ_i = a_i @ np.log(c_safe)
                r_i = lnQ_i - lnK[i]
                if abs(r_i) < tol:
                    continue

                # J_i = a_i @ diag(1/c) @ S  -> shape (R,)
                J_i = (a_i @ (inv_c[:, None] * S))
                grad_i = J_i * r_i

                step = learning_rate
                f_curr = 0.5 * (r_i * r_i)

                while True:
                    x_new = x - step * grad_i
                    c_new = c0 + S @ x_new
                    if np.all(c_new >= -1e-15):
                        c_new_safe = np.maximum(c_new, min_concentration)
                        lnQ_i_new = a_i @ np.log(c_new_safe)
                        r_i_new = lnQ_i_new - lnK[i]
                        f_new = 0.5 * (r_i_new * r_i_new)
                        if f_new <= f_curr or step < 1e-12:
                            x = x_new
                            any_update = True
                            break
                    step *= backtrack_beta

            # Full residual check for convergence
            c_full = c0 + S @ x
            c_full_safe = np.maximum(c_full, min_concentration)
            full_residual = (A @ np.log(c_full_safe)) - lnK
            if np.linalg.norm(full_residual, ord=2) < tol:
                break
            if not any_update:
                break

        c_final = c0 + S @ x
        c_final = np.maximum(c_final, 0.0)

        self.x_solution = x
        self.fitted = True
        return c_final.tolist()

    def _calculate_by_newton(self,
                              max_iter: int = 200,
                              learning_rate: float = 1.0,
                              tol: float = 1e-10,
                              backtrack_beta: float = 0.5,
                              min_concentration: float = 1e-12):
        # Stoichiometry and mappings
        N = self.env.stoichiometric_coefficient_array
        S = N.T
        R = N.shape[0]
        C = N.shape[1]

        c0 = np.array(self.env.concentrations, dtype=float)

        # Mass-action exponents, exclude s/l phases
        A = -N.astype(float)
        phase_include_mask = np.ones(C, dtype=bool)
        for j, compound in enumerate(self.env.compounds):
            ph = compound.phase(self.env.T)
            if ph in ("s", "l"):
                phase_include_mask[j] = False
        for j in range(C):
            if not phase_include_mask[j]:
                A[:, j] = 0.0

        K_vec = np.array([max(rxn.K, 1e-300) for rxn in self.env.reactions], dtype=float)
        lnK = np.log(K_vec)

        x = np.zeros(R, dtype=float)

        for _ in range(max_iter):
            c = c0 + S @ x
            c_safe = np.maximum(c, min_concentration)
            lnQ = A @ np.log(c_safe)
            r = lnQ - lnK
            if np.linalg.norm(r, ord=2) < tol:
                break

            inv_c = 1.0 / c_safe
            J = A @ (inv_c[:, None] * S)

            # Solve J * dx = r, then x <- x - alpha * dx
            try:
                dx, *_ = np.linalg.lstsq(J, r, rcond=None)
            except Exception:
                dx = np.linalg.pinv(J) @ r

            step = learning_rate
            f_curr = 0.5 * np.dot(r, r)
            while True:
                x_new = x - step * dx
                c_new = c0 + S @ x_new
                if np.all(c_new >= -1e-15):
                    c_new_safe = np.maximum(c_new, min_concentration)
                    lnQ_new = A @ np.log(c_new_safe)
                    r_new = lnQ_new - lnK
                    f_new = 0.5 * np.dot(r_new, r_new)
                    if f_new <= f_curr or step < 1e-12:
                        x = x_new
                        break
                step *= backtrack_beta

        c_final = c0 + S @ x
        c_final = np.maximum(c_final, 0.0)

        self.x_solution = x
        self.fitted = True
        return c_final.tolist()