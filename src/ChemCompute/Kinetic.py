import warnings

from ._general import Enviroment
from ._kinetics import integrate_kinetics


class KineticalCalculator:
    """
    Deprecated wrapper around ``Enviroment.kinetics(...)``.

    Prefer calling ``env.kinetics(...)`` directly on an ``Enviroment`` instance.
    """

    def __init__(self, accuracy=1e-3):
        warnings.warn(
            "KineticalCalculator is deprecated; use env.kinetics(...) instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        self.accuracy = accuracy
        self.fitted = False
        self.enviroment = None

    def fit(self, enviroment):
        if isinstance(enviroment, Enviroment):
            self.enviroment = enviroment
        else:
            raise ValueError("The input should be an instance of Enviroment class")
        self.rate_constants = enviroment.rate_constants
        self.reactions_by_index = enviroment.reaction_by_index
        self.stoichiometric_coefficient_by_reaction = enviroment.stoichiometric_coefficient_by_reaction
        self.rate_dependency_by_reaction = enviroment.rate_dependency_by_reaction
        self.number_of_reactions = len(enviroment)
        self.concentrations = [
            compound["concentration"] for compound in enviroment.compounds_concentration
        ]
        self.fitted = True

    def calculate(self, time, checkpoint_time=None, plot=False, directory="./plot.png", colors=None):
        if not self.fitted:
            raise NameError("You should fit the model to an enviromt object before calculation")
        if checkpoint_time is None:
            checkpoint_time = []
        return self.enviroment.kinetics(
            time=time,
            checkpoint_time=checkpoint_time,
            plot=plot,
            directory=directory,
            colors=colors,
            accuracy=self.accuracy,
        )

    def fit_calculate(
        self, enviroment, time, checkpoint_time=None, plot=False, directory="./plot.png", colors=None
    ):
        if checkpoint_time is None:
            checkpoint_time = []
        self.fit(enviroment)
        return self.calculate(time, checkpoint_time, plot, directory, colors)

    def calculate_responsively(self, checkpoint_time=None, animation_update_interval=0.1, colors=None):
        import matplotlib
        import matplotlib.pyplot as plt
        import numpy as np
        import random
        from itertools import count
        from matplotlib.animation import FuncAnimation

        if checkpoint_time is None:
            checkpoint_time = []

        if not self.fitted:
            raise NameError("You must fit the model to an Enviroment before calculation.")

        matplotlib.use("TkAgg", force=True)

        plt.figure()
        plot_colors = []
        checkpoints = []
        num_compounds = len(self.enviroment.compounds)

        if colors is not None:
            if len(colors) != num_compounds:
                raise ValueError(
                    f"Number of colors ({len(colors)}) must equal number of compounds ({num_compounds})"
                )
            plot_colors = colors
        else:
            for _ in self.enviroment.compounds:
                plot_colors.append(
                    (
                        random.randint(0, 95) / 100,
                        random.randint(0, 95) / 100,
                        random.randint(0, 95) / 100,
                    )
                )

        plt.xlabel("time")
        plt.ylabel("concentration")

        concentrations = self.enviroment.concentrations_array
        rate_dependencies = self.enviroment.rate_dependency_array
        stoichiometric_coefficient = self.enviroment.stoichiometric_coefficient_array
        rate_constants = self.enviroment.rate_constants_array
        time_interval = self.accuracy
        eps = 1e-300

        def calculate_rf():
            log_c = np.log(concentrations + eps)
            log_prod = rate_dependencies[:, 0, :] @ log_c
            return np.exp(log_prod) * rate_constants[:, 0] * time_interval

        def calculate_rb():
            log_c = np.log(concentrations + eps)
            log_prod = rate_dependencies[:, 1, :] @ log_c
            return np.exp(log_prod) * rate_constants[:, 1] * time_interval

        def calculate_concentration_change():
            rate = -(calculate_rf() - calculate_rb())
            return (stoichiometric_coefficient.T @ rate).reshape(-1)

        step_counter = count()

        def animate(_i):
            nonlocal concentrations
            t = self.accuracy * next(step_counter)

            new_concentrations = np.add(concentrations, calculate_concentration_change())
            new_concentrations[new_concentrations < 0] = 0

            for k in range(len(self.concentrations)):
                plt.plot(
                    [t, t - self.accuracy],
                    [new_concentrations[k], concentrations[k]],
                    color=plot_colors[k],
                )

            for checkpoint_t in checkpoint_time:
                if t <= checkpoint_t < t + self.accuracy:
                    checkpoints.append(new_concentrations.copy())
            concentrations = new_concentrations

        ani = FuncAnimation(
            plt.gcf(), animate, interval=animation_update_interval, cache_frame_data=False
        )

        for k in range(len(self.concentrations)):
            plt.plot(
                [0, 0],
                [0, 0],
                color=plot_colors[k],
                label=self.enviroment.compounds[k].unicode_formula,
            )

        plt.legend()
        plt.show(block=False)
        print("Type 'exit' to close / 'stop' to pause / 'resume' to continue:")
        while True:
            cmd = input().strip().lower()
            if cmd == "exit":
                ani.pause()
                plt.close()
                break
            if cmd == "stop":
                ani.pause()
            elif cmd == "resume":
                ani.resume()
            else:
                print("invalid_input")

        checkpoints.append(concentrations)
        return checkpoints
