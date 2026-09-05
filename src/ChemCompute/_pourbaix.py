"""Pourbaix diagram scanning via fixed pH and electrode potential."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

from ._half_reaction import BoundaryLine, compute_pH
from ._titration import _find_h_plus_index


def _pH_grid(pH_min: float, pH_max: float, steps: int) -> np.ndarray:
    if steps < 2:
        raise ValueError("pH grid steps must be at least 2.")
    if pH_min > pH_max:
        raise ValueError("Require pH_min <= pH_max.")
    return np.linspace(pH_min, pH_max, int(steps))


def _eh_grid(eh_min: float, eh_max: float, steps: int) -> np.ndarray:
    if steps < 2:
        raise ValueError("Eh grid steps must be at least 2.")
    if eh_min > eh_max:
        raise ValueError("Require Eh_min <= Eh_max.")
    return np.linspace(eh_min, eh_max, int(steps))


def _dominant_label(env, concentrations, track_species: Sequence[str]) -> str:
    label_map = {compound.formula: j for j, compound in enumerate(env.compounds)}
    best_label = track_species[0]
    best_value = -1.0
    for label in track_species:
        idx = label_map.get(label)
        if idx is None:
            continue
        value = float(concentrations[idx])
        if value > best_value:
            best_value = value
            best_label = label
    return best_label


@dataclass
class PourbaixResult:
    grid_pH: np.ndarray
    grid_Eh: np.ndarray
    grid_dominant: np.ndarray
    boundary_lines: list[BoundaryLine] = field(default_factory=list)
    speciation: Optional[np.ndarray] = None
    track_species: list[str] = field(default_factory=list)

    def matrix(self) -> np.ndarray:
        """Return dominant-species index grid (rows pH, cols Eh)."""
        return self.grid_dominant

    def plot_predominance(
        self,
        *,
        ax=None,
        show: bool = True,
        save: Optional[str] = None,
        labels: Optional[dict[int, str]] = None,
    ):
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(8, 6))
        unique = sorted(set(self.grid_dominant.ravel().tolist()))
        label_map = labels or {i: self.track_species[i] for i in range(len(self.track_species))}
        cmap = plt.get_cmap("tab20", max(len(unique), 1))
        image = ax.imshow(
            self.grid_dominant,
            origin="lower",
            aspect="auto",
            extent=[self.grid_Eh[0], self.grid_Eh[-1], self.grid_pH[0], self.grid_pH[-1]],
            cmap=cmap,
        )
        ax.set_xlabel("Eh (V vs SHE)")
        ax.set_ylabel("pH")
        ax.set_title("Pourbaix predominance")
        if save:
            plt.savefig(save, bbox_inches="tight")
        if show and not save:
            plt.show()
        return ax, image

    def plot_boundaries(
        self,
        *,
        ax=None,
        show: bool = True,
        save: Optional[str] = None,
        include_water_lines: bool = True,
    ):
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(8, 6))
        for line in self.boundary_lines:
            ax.plot(line.Eh_values, line.pH_values, label=line.label)
        if include_water_lines:
            pH = self.grid_pH
            eh_h2 = -0.05916 * pH
            eh_o2 = 1.229 - 0.05916 * pH
            ax.plot(eh_h2, pH, "--", color="gray", label="H+/H2")
            ax.plot(eh_o2, pH, "--", color="black", label="O2/H2O")
        ax.set_xlabel("Eh (V vs SHE)")
        ax.set_ylabel("pH")
        ax.set_title("Pourbaix boundaries")
        ax.legend(loc="best", fontsize=8)
        if save:
            plt.savefig(save, bbox_inches="tight")
        if show and not save:
            plt.show()
        return ax


@dataclass
class Pourbaix:
    """
    Pourbaix diagram scanner over pH and Eh using environment equilibrium.

    Each grid point fixes ``buffer=["H+"]``, sets ``electrode_Eh``, and calls
    ``equilibrium()`` on a copy of the base environment.
    """

    environment: object
    track_species: Sequence[str]
    pH_min: float = 0.0
    pH_max: float = 14.0
    pH_steps: int = 40
    Eh_min: float = -1.0
    Eh_max: float = 1.5
    Eh_steps: int = 40
    equilibrium_method: str = "newton"
    equilibrium_tol: float = 1e-8

    def run(self) -> PourbaixResult:
        if not self.track_species:
            raise ValueError("track_species must contain at least one formula label.")
        if _find_h_plus_index(self.environment) is None:
            raise ValueError("Pourbaix requires H+ in the environment.")

        pH_values = _pH_grid(self.pH_min, self.pH_max, self.pH_steps)
        eh_values = _eh_grid(self.Eh_min, self.Eh_max, self.Eh_steps)
        label_to_index = {label: i for i, label in enumerate(self.track_species)}
        dominant = np.zeros((len(pH_values), len(eh_values)), dtype=int)

        boundary_lines = []
        for hr in getattr(self.environment, "half_reactions", None) or []:
            boundary_lines.append(hr.boundary_line(pH_values))

        for i, pH in enumerate(pH_values):
            h_conc = 10.0 ** (-float(pH))
            for j, eh in enumerate(eh_values):
                env_copy = self.environment.copy()
                env_copy.set_buffer(["H+"])
                overrides = {"H+": h_conc}
                env_copy._apply_concentration_overrides(overrides)
                env_copy.set_electrode_potential(Eh=float(eh))
                try:
                    concentrations = env_copy.equilibrium(
                        method=self.equilibrium_method,
                        tol=self.equilibrium_tol,
                    )
                except Exception:
                    dominant[i, j] = label_to_index[self.track_species[0]]
                    continue
                label = _dominant_label(env_copy, concentrations, self.track_species)
                dominant[i, j] = label_to_index.get(label, 0)

        return PourbaixResult(
            grid_pH=pH_values,
            grid_Eh=eh_values,
            grid_dominant=dominant,
            boundary_lines=boundary_lines,
            track_species=list(self.track_species),
        )
