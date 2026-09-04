import warnings

from ._general import Enviroment


class EquilibriumCalculator:
    def __init__(self, method_of_calculation: str = "bgd"):
        warnings.warn(
            "EquilibriumCalculator is deprecated; use env.equilibrium(...) instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        self.method_of_calculation = method_of_calculation
        self.fitted = False
        self.env = None
        self.x_solution = None

    def _generate_concentration_equations(self):
        concentration_eq = [str(value) for value in self.env.concentrations]

        for r_index, reaction in enumerate(self.env.reactions, start=1):
            for idx, compound in enumerate(self.env.compounds):
                coeff = 0
                for reactant in reaction.reactants:
                    if reactant["compound"] == compound:
                        coeff += reactant["stoichiometric_coefficient"]
                        break
                for product in reaction.products:
                    if product["compound"] == compound:
                        coeff -= product["stoichiometric_coefficient"]
                        break
                if coeff != 0:
                    concentration_eq[idx] += f" + ({coeff}x{r_index})"
        return concentration_eq

    def fit(self, env: Enviroment):
        if isinstance(env, Enviroment):
            self.env = env
        else:
            raise ValueError("The input should be an instance of Enviroment class")
        self.concentration_equation = self._generate_concentration_equations()
        self.fitted = True

    def calculate(
        self,
        max_iter=None,
        learning_rate=None,
        tol=None,
        backtrack_beta: float = 0.5,
        min_concentration: float = 1e-12,
        loss: str = "log_quotient",
        concentration_error_limit=None,
        huber_delta: float = 1.0,
    ):
        if not self.fitted:
            raise ValueError("Environment not fitted")

        result = self.env.equilibrium(
            method=self.method_of_calculation,
            loss=loss,
            max_iter=max_iter,
            learning_rate=learning_rate,
            tol=tol,
            backtrack_beta=backtrack_beta,
            min_concentration=min_concentration,
            concentration_error_limit=concentration_error_limit,
            huber_delta=huber_delta,
        )
        self.x_solution = getattr(self.env, "_equilibrium_x_solution", None)
        return result

    def fit_calculate(
        self,
        env: Enviroment,
        max_iter=None,
        learning_rate=None,
        tol=None,
        backtrack_beta: float = 0.5,
        min_concentration: float = 1e-12,
        loss: str = "log_quotient",
        concentration_error_limit=None,
        huber_delta: float = 1.0,
    ):
        self.fit(env)
        return self.calculate(
            max_iter=max_iter,
            learning_rate=learning_rate,
            tol=tol,
            backtrack_beta=backtrack_beta,
            min_concentration=min_concentration,
            loss=loss,
            concentration_error_limit=concentration_error_limit,
            huber_delta=huber_delta,
        )
