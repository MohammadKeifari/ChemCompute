"""Titration curves via volume-aware environment mixing."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Optional, Sequence

import numpy as np


def _find_h_plus_index(env):
    for j, compound in enumerate(env.compounds):
        if compound.formula == "H+":
            return j
    return None


def _compute_pH(env, concentrations: np.ndarray) -> float:
    h_index = _find_h_plus_index(env)
    if h_index is None:
        return float("nan")
    h_conc = concentrations[h_index]
    if getattr(env, "activity_model", None) is not None:
        gammas = env.activity_model.gamma_array(env, concentrations, env.T)
        h_conc = gammas[h_index] * h_conc
    return -math.log10(max(h_conc, 1e-300))


def _detect_equivalence_points(pH_values: Sequence[float]) -> list[int]:
    """Indices of largest |dpH/dV| changes."""
    if len(pH_values) < 3:
        return []
    dpH = np.abs(np.diff(pH_values))
    if dpH.max() <= 0:
        return []
    threshold = 0.5 * dpH.max()
    return [i + 1 for i, val in enumerate(dpH) if val >= threshold]


def _titrant_slug(titrant, volume_added: float):
    """Build a titrant aliquot at the given volume from titrant stock concentrations."""
    from ._general import Enviroment

    buffer_spec = getattr(titrant, "_buffer_spec", None)
    return Enviroment.from_compounds(
        titrant.concentrations_dict,
        T=titrant.T,
        volume=float(volume_added),
        adjust_thermodynamics=titrant.adjust_thermodynamics,
        activity_model=titrant.activity_model,
        buffer=buffer_spec if buffer_spec else None,
    )


def mix_sample_with_titrant(sample, titrant, titrant_volume: float):
    """Combine sample with a titrant aliquot without mutating the originals."""
    from ._general import Enviroment

    sample_copy = sample.copy()
    if titrant_volume <= 0:
        return sample_copy
    if titrant.volume <= 0:
        raise ValueError("Titrant environment volume must be positive.")
    slug = _titrant_slug(titrant, titrant_volume)
    return Enviroment.combine((1.0, sample_copy), (1.0, slug))


def _volume_steps(
    volume_min: float,
    volume_max: float,
    steps: int,
    volumes: Optional[Sequence[float]] = None,
) -> np.ndarray:
    if volumes is not None:
        arr = np.asarray(volumes, dtype=float)
        if arr.ndim != 1 or len(arr) == 0:
            raise ValueError("volumes must be a non-empty 1D sequence.")
        if np.any(arr < 0):
            raise ValueError("Titrant volumes must be non-negative.")
        return arr
    if steps < 2:
        raise ValueError("steps must be at least 2 when volumes is not provided.")
    if volume_min < 0 or volume_max < volume_min:
        raise ValueError("Require 0 <= volume_min <= volume_max.")
    return np.linspace(volume_min, volume_max, int(steps))


@dataclass
class TitrationResult:
    """
    Titration curve: equilibrium concentrations vs titrant volume added.

    Attributes
    ----------
    titrant_volumes : list[float]
        Titrant volume added at each step (litres).
    total_volumes : list[float]
        Total solution volume after each addition (litres).
    compound_labels : list[str]
        Species labels aligned with concentration columns.
    concentration_matrix : np.ndarray
        Shape ``(n_steps, n_compounds)`` — rows indexed by ``titrant_volumes``.
    pH : list[float]
        pH at each step (NaN if H+ is absent).
    speciation : list[dict[str, float]]
        Per-step ``{formula: concentration}`` maps.
    equilibrium_results : list
        Full :class:`EquilibriumResult` objects per step.
    equivalence_hints : list[int]
        Step indices with large |dpH/dV| (informational).
    """

    titrant_volumes: list[float] = field(default_factory=list)
    total_volumes: list[float] = field(default_factory=list)
    compound_labels: list[str] = field(default_factory=list)
    concentration_matrix: np.ndarray = field(default_factory=lambda: np.empty((0, 0)))
    pH: list[float] = field(default_factory=list)
    speciation: list[dict[str, float]] = field(default_factory=list)
    equilibrium_results: list[Any] = field(default_factory=list)
    equivalence_hints: list[int] = field(default_factory=list)

    def matrix(self) -> np.ndarray:
        """
        Concentration matrix with shape ``(n_steps, n_compounds)``.

        Rows correspond to ``titrant_volumes``; columns to ``compound_labels``.
        """
        return np.array(self.concentration_matrix, copy=True)

    def species(self, formula: str) -> np.ndarray:
        """Concentration series for one species vs titrant volume added."""
        if formula not in self.compound_labels:
            raise KeyError(f"Species {formula!r} not in titration result.")
        index = self.compound_labels.index(formula)
        return self.concentration_matrix[:, index]

    def plot(
        self,
        species: Optional[Sequence[str]] = None,
        *,
        plot: bool | str = False,
        directory: str = "./titration.png",
        colors: Optional[Sequence[str]] = None,
    ):
        """
        Plot species concentrations vs titrant volume added.

        Parameters
        ----------
        species : sequence of str, optional
            Subset of ``compound_labels`` to plot. Default: all species.
        plot : False, ``"interactive"``, or ``"save"``
            Matplotlib output mode (same convention as ``env.kinetics()``).
        directory : str
            Output path when ``plot="save"``.
        colors : list, optional
            One color per plotted species.
        """
        if plot not in (False, "save", "interactive"):
            raise ValueError("`plot` is not one of [False, 'save', 'interactive'].")

        labels = list(species) if species is not None else list(self.compound_labels)
        if not labels:
            raise ValueError("No species to plot.")

        import matplotlib

        if plot == "interactive":
            matplotlib.use("TkAgg", force=True)
        elif plot == "save":
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        x = np.asarray(self.titrant_volumes, dtype=float)
        if colors is not None and len(colors) != len(labels):
            raise ValueError("Number of colors must equal number of species plotted.")

        for i, formula in enumerate(labels):
            kwargs = {"label": formula}
            if colors is not None:
                kwargs["color"] = colors[i]
            plt.plot(x, self.species(formula), **kwargs)

        plt.xlabel("Titrant volume added (L)")
        plt.ylabel("Concentration (mol/L)")
        plt.legend()

        if plot == "interactive":
            plt.show(block=False)
        elif plot == "save":
            plt.savefig(directory)
            plt.close("all")

    def plot_pH(
        self,
        *,
        plot: bool | str = False,
        directory: str = "./titration_pH.png",
        color: str = "#26547c",
    ):
        """Plot pH vs titrant volume added."""
        if plot not in (False, "save", "interactive"):
            raise ValueError("`plot` is not one of [False, 'save', 'interactive'].")
        if not self.pH or all(math.isnan(v) for v in self.pH):
            raise ValueError("pH data is not available for this titration.")

        import matplotlib

        if plot == "interactive":
            matplotlib.use("TkAgg", force=True)
        elif plot == "save":
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.plot(self.titrant_volumes, self.pH, color=color)
        plt.xlabel("Titrant volume added (L)")
        plt.ylabel("pH")
        plt.ylim(0, 14)

        if plot == "interactive":
            plt.show(block=False)
        elif plot == "save":
            plt.savefig(directory)
            plt.close("all")


class Titration:
    """
    Run a titration by mixing a sample environment with a titrant environment.

    At each step, ``titrant_volume`` litres of titrant (at its stock concentrations)
    are combined with the sample using volume-weighted mixing, then equilibrium
    is solved on the mixture.

    Parameters
    ----------
    sample : Enviroment
        Analyte / initial solution. Its ``volume`` sets the starting sample size.
    titrant : Enviroment
        Titrant stock (typically ``from_compounds``). Concentrations define titrant
        strength; ``volume`` defines the reference stock size (only concentrations
        matter for each aliquot).
    volume_min, volume_max, steps : float, float, int
        Titrant volume added at first step, last step, and number of steps
        (inclusive end points via ``numpy.linspace``).
    volumes : sequence of float, optional
        Explicit titrant volumes to use instead of ``volume_min``/``volume_max``/``steps``.
    """

    def __init__(
        self,
        sample,
        titrant,
        *,
        volume_min: float = 0.0,
        volume_max: float = 0.05,
        steps: int = 50,
        volumes: Optional[Sequence[float]] = None,
    ):
        from ._general import Enviroment

        if not isinstance(sample, Enviroment):
            raise TypeError("sample must be an Enviroment instance.")
        if not isinstance(titrant, Enviroment):
            raise TypeError("titrant must be an Enviroment instance.")
        if sample.volume <= 0:
            raise ValueError("Sample environment volume must be positive.")
        if titrant.volume <= 0:
            raise ValueError("Titrant environment volume must be positive.")
        if sample.T != titrant.T:
            raise ValueError("Sample and titrant must have the same temperature.")
        if sample.adjust_thermodynamics != titrant.adjust_thermodynamics:
            raise ValueError(
                "Sample and titrant must have the same adjust_thermodynamics setting."
            )

        self.sample = sample
        self.titrant = titrant
        self.volume_min = float(volume_min)
        self.volume_max = float(volume_max)
        self.steps = int(steps)
        self.volumes = volumes

    def run(self, **equilibrium_kwargs) -> TitrationResult:
        """Execute the titration and return a :class:`TitrationResult`."""
        volume_steps = _volume_steps(
            self.volume_min,
            self.volume_max,
            self.steps,
            self.volumes,
        )

        result = TitrationResult()
        step_specs = []
        step_meta = []

        for titrant_volume in volume_steps:
            mixed = mix_sample_with_titrant(self.sample, self.titrant, float(titrant_volume))
            eq_result = mixed.equilibrium(return_details=True, **equilibrium_kwargs)
            conc = np.array(eq_result.concentrations, dtype=float)
            pH = _compute_pH(mixed, conc)
            spec = dict(zip(mixed.compound_labels, eq_result.concentrations))

            step_specs.append(spec)
            step_meta.append(
                {
                    "titrant_volume": float(titrant_volume),
                    "total_volume": float(mixed.volume),
                    "pH": pH,
                    "eq_result": eq_result,
                }
            )

        label_order = []
        seen = set()
        for spec in step_specs:
            for label in spec:
                if label not in seen:
                    seen.add(label)
                    label_order.append(label)

        rows = [[spec.get(label, 0.0) for label in label_order] for spec in step_specs]
        for meta, spec in zip(step_meta, step_specs):
            result.titrant_volumes.append(meta["titrant_volume"])
            result.total_volumes.append(meta["total_volume"])
            result.pH.append(meta["pH"])
            result.speciation.append(spec)
            result.equilibrium_results.append(meta["eq_result"])

        result.compound_labels = label_order
        if rows:
            result.concentration_matrix = np.asarray(rows, dtype=float)
        result.equivalence_hints = _detect_equivalence_points(result.pH)
        return result
