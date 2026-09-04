import random

import matplotlib
import numpy as np


def integrate_kinetics(
    env,
    time,
    accuracy=1e-3,
    checkpoint_time=None,
    plot=False,
    directory="./plot.png",
    colors=None,
):
    """
    Numerically integrate reaction kinetics for an environment.

    Parameters
    ----------
    env : Enviroment
        Reaction environment.
    time : float
        Total simulation time.
    accuracy : float
        Integration time step.
    checkpoint_time : list[float], optional
        Times at which to record concentrations.
    plot : bool or str
        ``False``, ``"interactive"``, or ``"save"``.
    directory : str
        Output path when ``plot="save"``.
    colors : list, optional
        Plot colors, one per compound.

    Returns
    -------
    list
        Checkpoint concentration snapshots.
    """
    if checkpoint_time is None:
        checkpoint_time = []

    if plot not in (False, "save", "interactive"):
        raise ValueError("`plot` is not one of [False, 'save', 'interactive'].")

    if plot == "interactive":
        matplotlib.use("TkAgg", force=True)
    elif plot == "save":
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_colors = []
    checkpoints = []
    if plot:
        plt.figure()
        num_compounds = len(env.compounds)

        if colors is not None:
            if len(colors) != num_compounds:
                raise ValueError(
                    f"Number of colors ({len(colors)}) must equal number of compounds ({num_compounds})"
                )
            plot_colors = colors
        else:
            for _ in env.compounds:
                plot_colors.append(
                    (
                        random.randint(0, 95) / 100,
                        random.randint(0, 95) / 100,
                        random.randint(0, 95) / 100,
                    )
                )

        plt.xlabel("time")
        plt.ylabel("concentration")

    concentrations = env.concentrations_array.copy()
    rate_dependencies = env.rate_dependency_array
    stoichiometric_coefficient = env.stoichiometric_coefficient_array
    rate_constants = env.rate_constants_array
    time_interval = accuracy
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

    t = 0.0
    for _ in range(int(time / accuracy + 1)):
        new_concentrations = np.add(concentrations, calculate_concentration_change())
        new_concentrations[new_concentrations < 0] = 0
        if plot:
            for k in range(num_compounds):
                plt.plot(
                    [t, t - accuracy],
                    [new_concentrations[k], concentrations[k]],
                    color=plot_colors[k],
                )
        for checkpoint_t in checkpoint_time:
            if t <= checkpoint_t < t + accuracy:
                checkpoints.append(new_concentrations.copy())
        concentrations = new_concentrations
        t += accuracy

    if plot == "interactive":
        for k in range(num_compounds):
            plt.plot(
                [0, 0],
                [0, 0],
                color=plot_colors[k],
                label=env.compounds[k].unicode_formula,
            )
        plt.legend()
        plt.show(block=False)

        print("Type 'exit' to close the plot:")
        while True:
            cmd = input().strip().lower()
            if cmd == "exit":
                plt.close()
                break
            print("Invalid input.")
    elif plot == "save":
        for k in range(num_compounds):
            plt.plot(
                [0, 0],
                [0, 0],
                color=plot_colors[k],
                label=env.compounds[k].unicode_formula,
            )
        plt.legend()
        plt.savefig(directory)
        plt.close("all")
        del plt

    checkpoints.append(concentrations)
    return checkpoints
