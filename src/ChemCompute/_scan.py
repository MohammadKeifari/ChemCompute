"""Parameter scans: titration curves and Pourbaix diagrams."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Optional, Sequence

import numpy as np

from ._buffer import buffer_diagnostics


FARADAY = 96485.33212
R_GAS = 8.3145


@dataclass
class ScanResult:
    """Results from a parameter scan."""

    axis: str
    x_values: list[float] = field(default_factory=list)
    pH: list[float] = field(default_factory=list)
    concentrations: list[list[float]] = field(default_factory=list)
    speciation: list[dict[str, float]] = field(default_factory=list)
    equilibrium_results: list[Any] = field(default_factory=list)
    equivalence_hints: list[int] = field(default_factory=list)
    dominant_species: list[str] = field(default_factory=list)
    grid_pH: Optional[np.ndarray] = None
    grid_Eh: Optional[np.ndarray] = None
    grid_dominant: Optional[np.ndarray] = None
    absorbance: Optional[np.ndarray] = None
    wavelengths: Optional[np.ndarray] = None


def _find_h_plus_index(env) -> Optional[int]:
    for j, compound in enumerate(env.compounds):
        if compound.formula in ("H+", "H+"):
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


def _apply_titrant_dilution(
    env,
    sample_volume: float,
    titrant_volume: float,
    titrant_concentration: float,
    titrant_formula: str,
) -> None:
    """Dilute sample and add titrant moles."""
    total_volume = sample_volume + titrant_volume
    scale = sample_volume / total_volume
    new_conc = [c * scale for c in env.concentrations]

    if titrant_formula in env.concentrations_dict:
        idx = env.compound_labels.index(titrant_formula)
        added = titrant_concentration * titrant_volume / total_volume
        new_conc[idx] += added
    else:
        raise ValueError(
            f"Titrant species {titrant_formula!r} not in environment. "
            "Add a dissociation reaction or include the titrant species."
        )
    env.concentrations = new_conc


def _detect_equivalence_points(pH_values: Sequence[float]) -> list[int]:
    """Indices of largest |dpH/dV| changes."""
    if len(pH_values) < 3:
        return []
    dpH = np.abs(np.diff(pH_values))
    if dpH.max() <= 0:
        return []
    threshold = 0.5 * dpH.max()
    hints = [i + 1 for i, val in enumerate(dpH) if val >= threshold]
    return hints


def _apply_fixed_pH(env, pH: float) -> None:
    h_index = _find_h_plus_index(env)
    if h_index is None:
        raise ValueError("Pourbaix scan requires H+ in the environment.")
    h_conc = 10.0 ** (-pH)
    conc = list(env.concentrations)
    conc[h_index] = h_conc
    env.concentrations = conc


def _apply_redox_potential(env, redox_couples: list[dict], Eh: float, pH: float) -> None:
    """Adjust equilibrium K of redox reactions from Nernst relation (vs SHE)."""
    for couple in redox_couples:
        reaction_index = couple["reaction_index"]
        n_electrons = couple.get("n_electrons", 1)
        E0 = couple["E0"]
        rxn = env.reactions[reaction_index]
        # Nernst: E = E0 - (0.05916/n)*pH at standard; K = exp(nF(Eh - E)/RT)
        E_effective = E0 - (0.05916 / n_electrons) * pH
        delta = n_electrons * FARADAY * (Eh - E_effective) / (R_GAS * env.T)
        rxn.K = math.exp(delta)


class ParameterScan:
    """
    Orchestrate equilibrium solves over a scan axis.

    Supports titration (``axis='titrant_volume'``) and Pourbaix grids (``axis='grid'``).
    """

    def __init__(
        self,
        base_env,
        axis: str,
        *,
        titrant: Optional[dict] = None,
        sample_volume: float = 0.1,
        pH_range: Optional[tuple[float, float, float]] = None,
        Eh_range: Optional[tuple[float, float, float]] = None,
        redox_couples: Optional[list[dict]] = None,
        track_species: Optional[list[str]] = None,
        uvvis_wavelengths: Optional[Sequence[float]] = None,
        uvvis_path_length: float = 0.01,
    ):
        self.base_env = base_env
        self.axis = axis
        self.titrant = titrant or {}
        self.sample_volume = sample_volume
        self.pH_range = pH_range
        self.Eh_range = Eh_range
        self.redox_couples = redox_couples or []
        self.track_species = track_species
        self.uvvis_wavelengths = uvvis_wavelengths
        self.uvvis_path_length = uvvis_path_length

    def run_equilibrium(self, **equilibrium_kwargs) -> ScanResult:
        if self.axis == "titrant_volume":
            return self._run_titration(**equilibrium_kwargs)
        if self.axis == "grid":
            return self._run_pourbaix(**equilibrium_kwargs)
        raise ValueError("axis must be 'titrant_volume' or 'grid'")

    def _run_titration(self, **equilibrium_kwargs) -> ScanResult:
        titrant_formula = self.titrant.get("formula", "H+")
        titrant_conc = self.titrant.get("concentration", 0.1)
        volume_steps = self.titrant.get("volume_steps")
        if volume_steps is None:
            volume_steps = np.linspace(0.0, 0.05, 100)
        volume_steps = np.asarray(volume_steps, dtype=float)

        result = ScanResult(axis="titrant_volume")
        from ._uvvis import uvvis_spectrum

        absorbance_rows = []

        for vol in volume_steps:
            env = self.base_env.copy()
            _apply_titrant_dilution(
                env,
                self.sample_volume,
                float(vol),
                titrant_conc,
                titrant_formula,
            )
            eq_result = env.equilibrium(return_details=True, **equilibrium_kwargs)
            conc = np.array(eq_result.concentrations, dtype=float)
            pH = _compute_pH(env, conc)

            result.x_values.append(float(vol))
            result.pH.append(pH)
            result.concentrations.append(eq_result.concentrations)
            spec = dict(zip(env.compound_labels, eq_result.concentrations))
            result.speciation.append(spec)
            result.equilibrium_results.append(eq_result)

            if self.track_species:
                dominant = max(
                    ((s, spec.get(s, 0.0)) for s in self.track_species),
                    key=lambda item: item[1],
                )[0]
            else:
                dominant = max(spec.items(), key=lambda item: item[1])[0]
            result.dominant_species.append(dominant)

            if self.uvvis_wavelengths is not None or getattr(env, "spectra", None):
                wl = self.uvvis_wavelengths
                absorbance_rows.append(
                    uvvis_spectrum(
                        env,
                        eq_result.concentrations,
                        wl,
                        path_length=self.uvvis_path_length,
                    )
                )

        result.equivalence_hints = _detect_equivalence_points(result.pH)
        if absorbance_rows:
            result.wavelengths = np.asarray(
                self.uvvis_wavelengths
                or sorted(
                    {
                        wl
                        for spec in getattr(self.base_env, "spectra", {}).values()
                        for wl, _ in spec.points
                    }
                ),
                dtype=float,
            )
            result.absorbance = np.vstack(absorbance_rows)
        return result

    def _run_pourbaix(self, **equilibrium_kwargs) -> ScanResult:
        if self.pH_range is None or self.Eh_range is None:
            raise ValueError("Pourbaix scan requires pH_range and Eh_range.")
        pH_start, pH_stop, pH_step = self.pH_range
        Eh_start, Eh_stop, Eh_step = self.Eh_range
        pH_values = np.arange(pH_start, pH_stop + 0.5 * pH_step, pH_step)
        Eh_values = np.arange(Eh_start, Eh_stop + 0.5 * Eh_step, Eh_step)

        result = ScanResult(axis="grid")
        result.grid_pH = pH_values
        result.grid_Eh = Eh_values
        grid_dominant = np.empty((len(Eh_values), len(pH_values)), dtype=object)

        for i, Eh in enumerate(Eh_values):
            for j, pH in enumerate(pH_values):
                env = self.base_env.copy()
                # Restore redox K from reference copies
                for couple in self.redox_couples:
                    rxn = env.reactions[couple["reaction_index"]]
                    rxn.K = couple["K_ref"]
                _apply_fixed_pH(env, float(pH))
                if self.redox_couples:
                    _apply_redox_potential(env, self.redox_couples, float(Eh), float(pH))
                eq_result = env.equilibrium(return_details=True, **equilibrium_kwargs)
                spec = dict(zip(env.compound_labels, eq_result.concentrations))
                if self.track_species:
                    dominant = max(
                        ((s, spec.get(s, 0.0)) for s in self.track_species),
                        key=lambda item: item[1],
                    )[0]
                else:
                    dominant = max(spec.items(), key=lambda item: item[1])[0]
                grid_dominant[i, j] = dominant

        result.grid_dominant = grid_dominant
        return result


def prepare_redox_couple(env, reaction_index: int, E0: float, n_electrons: int = 1) -> dict:
    """Build a redox couple descriptor for Pourbaix scans."""
    return {
        "reaction_index": reaction_index,
        "E0": E0,
        "n_electrons": n_electrons,
        "K_ref": env.reactions[reaction_index].K,
    }
