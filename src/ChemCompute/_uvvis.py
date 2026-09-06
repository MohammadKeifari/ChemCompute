"""UV-Vis Beer-Lambert spectra from speciation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional, Sequence, Union

import numpy as np


@dataclass
class SpectrumSpec:
    """
    Piecewise molar absorptivity specification.

    Parameters
    ----------
    points : sequence of (wavelength, epsilon)
        Wavelength in metres; epsilon in M^-1 m^-1 (or consistent path-length units).
    extrapolate : str
        ``flat`` — linear between points; flat beyond endpoints.
        ``none`` — epsilon is zero except at exactly listed wavelengths.
    """

    points: Sequence[tuple[float, float]]
    extrapolate: str = "flat"

    def __post_init__(self):
        if self.extrapolate not in ("flat", "none"):
            raise ValueError("extrapolate must be 'flat' or 'none'")
        sorted_points = sorted(self.points, key=lambda item: item[0])
        self.points = list(sorted_points)

    def epsilon(self, wavelength: float) -> float:
        """Return molar absorptivity at a single wavelength."""
        if not self.points:
            return 0.0
        wl = float(wavelength)
        if self.extrapolate == "none":
            for point_wl, eps in self.points:
                if np.isclose(wl, point_wl):
                    return eps
            return 0.0

        wl_list = [p[0] for p in self.points]
        eps_list = [p[1] for p in self.points]
        if wl <= wl_list[0]:
            return eps_list[0]
        if wl >= wl_list[-1]:
            return eps_list[-1]
        for idx in range(len(wl_list) - 1):
            if wl_list[idx] <= wl <= wl_list[idx + 1]:
                t = (wl - wl_list[idx]) / (wl_list[idx + 1] - wl_list[idx])
                return eps_list[idx] + t * (eps_list[idx + 1] - eps_list[idx])
        return 0.0

    def epsilon_array(self, wavelengths: Iterable[float]) -> np.ndarray:
        return np.array([self.epsilon(w) for w in wavelengths], dtype=float)


def _compound_spectrum(compound) -> Optional[SpectrumSpec]:
    return getattr(compound, "spectrum", None)


def _resolve_spectrum(compound, spectra: Optional[dict] = None) -> Optional[SpectrumSpec]:
    formula = compound.formula if hasattr(compound, "formula") else str(compound)
    if spectra is not None:
        return spectra.get(formula)
    return _compound_spectrum(compound)


def uvvis_spectrum(
    env,
    concentrations=None,
    wavelengths=None,
    *,
    path_length: float = 0.01,
    spectra: Optional[dict] = None,
) -> np.ndarray:
    """
    Compute absorbance A(λ) = sum_i epsilon_i(λ) * c_i * l.

    Molar absorptivity is read from each :class:`Compound`'s ``spectrum``.
    Pass ``spectra={formula: SpectrumSpec, ...}`` to override for this call.

    Parameters
    ----------
    path_length : float
        Path length in metres (default 0.01 m = 1 cm).
    """
    if concentrations is None:
        conc = np.array(env.concentrations, dtype=float)
    else:
        conc = np.array(concentrations, dtype=float)

    if wavelengths is None:
        all_wl = set()
        if spectra is not None:
            specs = list(spectra.values())
        else:
            specs = [
                compound.spectrum
                for compound in env.compounds
                if getattr(compound, "spectrum", None) is not None
            ]
        for spec in specs:
            for wl, _ in spec.points:
                all_wl.add(wl)
        if not all_wl:
            raise ValueError("No wavelengths provided and no compound spectra registered.")
        wavelengths = sorted(all_wl)

    wavelengths = np.asarray(list(wavelengths), dtype=float)
    absorbance = np.zeros_like(wavelengths, dtype=float)

    for j, compound in enumerate(env.compounds):
        spec = _resolve_spectrum(compound, spectra)
        if spec is None:
            continue
        eps = spec.epsilon_array(wavelengths)
        absorbance += eps * conc[j] * path_length

    return absorbance
