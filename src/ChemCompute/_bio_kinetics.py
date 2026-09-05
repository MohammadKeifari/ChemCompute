"""Enzyme kinetics integration (Michaelis-Menten and inhibition)."""

from __future__ import annotations

import random
from typing import Optional

import matplotlib
import numpy as np


def _mm_rate(concentrations, env, rxn):
    """Michaelis-Menten rate v = Vmax[S]/(Km+[S])."""
    params = rxn.bio_params
    substrate = params["substrate"]
    s_idx = env.compound_labels.index(substrate)
    s_conc = max(concentrations[s_idx], 0.0)
    vmax = params["Vmax"]
    km = params["Km"]
    return vmax * s_conc / (km + s_conc + 1e-300)


def _mm_competitive_rate(concentrations, env, rxn):
    params = rxn.bio_params
    substrate = params["substrate"]
    inhibitor = params["inhibitor"]
    s_idx = env.compound_labels.index(substrate)
    i_idx = env.compound_labels.index(inhibitor)
    s_conc = max(concentrations[s_idx], 0.0)
    i_conc = max(concentrations[i_idx], 0.0)
    vmax = params["Vmax"]
    km = params["Km"]
    ki = params["Ki"]
    apparent_km = km * (1.0 + i_conc / ki)
    return vmax * s_conc / (apparent_km + s_conc + 1e-300)


def _mm_uncompetitive_rate(concentrations, env, rxn):
    params = rxn.bio_params
    substrate = params["substrate"]
    inhibitor = params["inhibitor"]
    s_idx = env.compound_labels.index(substrate)
    i_idx = env.compound_labels.index(inhibitor)
    s_conc = max(concentrations[s_idx], 0.0)
    i_conc = max(concentrations[i_idx], 0.0)
    vmax = params["Vmax"]
    km = params["Km"]
    ki = params["Ki"]
    apparent_vmax = vmax / (1.0 + i_conc / ki)
    apparent_km = km / (1.0 + i_conc / ki)
    return apparent_vmax * s_conc / (apparent_km + s_conc + 1e-300)


def _mm_noncompetitive_rate(concentrations, env, rxn):
    params = rxn.bio_params
    substrate = params["substrate"]
    inhibitor = params["inhibitor"]
    s_idx = env.compound_labels.index(substrate)
    i_idx = env.compound_labels.index(inhibitor)
    s_conc = max(concentrations[s_idx], 0.0)
    i_conc = max(concentrations[i_idx], 0.0)
    vmax = params["Vmax"]
    km = params["Km"]
    ki = params["Ki"]
    apparent_vmax = vmax / (1.0 + i_conc / ki)
    return apparent_vmax * s_conc / (km + s_conc + 1e-300)


def _mm_mixed_rate(concentrations, env, rxn):
    params = rxn.bio_params
    substrate = params["substrate"]
    inhibitor = params["inhibitor"]
    s_idx = env.compound_labels.index(substrate)
    i_idx = env.compound_labels.index(inhibitor)
    s_conc = max(concentrations[s_idx], 0.0)
    i_conc = max(concentrations[i_idx], 0.0)
    vmax = params["Vmax"]
    km = params["Km"]
    ki = params["Ki"]
    alpha = params.get("alpha", 1.0)
    alpha_prime = params.get("alpha_prime", 1.0)
    apparent_vmax = vmax / (1.0 + i_conc / (alpha_prime * ki))
    apparent_km = km * (1.0 + i_conc / (alpha * ki))
    return apparent_vmax * s_conc / (apparent_km + s_conc + 1e-300)


RATE_LAW_FUNCTIONS = {
    "michaelis_menten": _mm_rate,
    "mm_competitive": _mm_competitive_rate,
    "mm_uncompetitive": _mm_uncompetitive_rate,
    "mm_noncompetitive": _mm_noncompetitive_rate,
    "mm_mixed": _mm_mixed_rate,
}


def _bio_reaction_rate(concentrations, env, rxn, direction: str = "forward"):
    rate_law = getattr(rxn, "rate_law", "mass_action")
    if rate_law not in RATE_LAW_FUNCTIONS:
        raise ValueError(f"Unknown bio rate law: {rate_law}")
    rate = RATE_LAW_FUNCTIONS[rate_law](concentrations, env, rxn)
    return rate if direction == "forward" else 0.0


def integrate_bio_kinetics(
    env,
    time,
    accuracy=1e-3,
    checkpoint_time=None,
    plot=False,
    directory="./plot.png",
    colors=None,
):
    """
    Integrate enzyme kinetics using MM and inhibition rate laws.

    Bio reactions produce net rate on their product species via stoichiometry.
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
                raise ValueError("Number of colors must equal number of compounds.")
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
    stoichiometric_coefficient = env.stoichiometric_coefficient_array
    time_interval = accuracy

    def net_rate_vector():
        rates = np.zeros(len(env.compounds), dtype=float)
        for rxn_index, rxn in enumerate(env.reactions):
            rate_law = getattr(rxn, "rate_law", "mass_action")
            if rate_law in RATE_LAW_FUNCTIONS:
                v = _bio_reaction_rate(concentrations, env, rxn, "forward")
            else:
                eps = 1e-300
                log_c = np.log(concentrations + eps)
                rate_dependencies = env.rate_dependency_array
                rate_constants = env.rate_constants_array
                log_prod_f = rate_dependencies[rxn_index, 0, :] @ log_c
                log_prod_b = rate_dependencies[rxn_index, 1, :] @ log_c
                rf = np.exp(log_prod_f) * rate_constants[rxn_index, 0]
                rb = np.exp(log_prod_b) * rate_constants[rxn_index, 1]
                v = rf - rb
            rates += stoichiometric_coefficient[rxn_index, :] * (-v)
        return rates

    t = 0.0
    for _ in range(int(time / accuracy + 1)):
        dc = net_rate_vector() * time_interval
        new_concentrations = concentrations + dc
        new_concentrations[new_concentrations < 0] = 0.0
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

    checkpoints.append(concentrations)
    return checkpoints


def uses_bio_kinetics(env) -> bool:
    return any(getattr(rxn, "rate_law", "mass_action") in RATE_LAW_FUNCTIONS for rxn in env.reactions)
