"""Shared constructors for the compound library."""

from __future__ import annotations

import math

from .._general import Compound
from .._uvvis import SpectrumSpec

T_REF = 298.0
_LN2_4 = 4.0 * math.log(2.0)


def library(fn):
    """Cache a library object so every call returns the same instance."""
    cached = []

    def getter():
        if not cached:
            cached.append(fn())
        return cached[0]

    getter.__name__ = fn.__name__
    getter.__qualname__ = getattr(fn, "__qualname__", fn.__name__)
    getter.__doc__ = fn.__doc__
    getter.__module__ = fn.__module__
    return getter


def phase_at(phase, temperature=T_REF):
    return [{"phase": phase, "temperature": temperature}]


def vis_nm(*peaks, extrapolate="flat"):
    """Polyline spectrum from ``(wavelength_nm, epsilon_M_minus_1_cm_minus_1)``.

    Points are stored as metres and M^-1 m^-1 so
    ``uvvis_spectrum(..., path_length=0.01)`` matches a 1 cm cuvette.
    ``extrapolate="flat"`` linearly connects neighbouring points. Put ε = 0
    at the wings so the curve falls to zero off-band.
    """
    if not peaks:
        raise ValueError("need at least one (nm, epsilon) point")
    return SpectrumSpec(
        points=[(float(nm) * 1e-9, float(eps) * 100.0) for nm, eps in peaks],
        extrapolate=extrapolate,
    )


def vis_bands(*bands, fwhm_cm=3500.0, step_nm=10.0):
    """Connected envelope from one or more Gaussian bands in wavenumber.

    Each band is ``(lambda_nm, epsilon_max)`` or
    ``(lambda_nm, epsilon_max, fwhm_cm_minus_1)``. Overlapping bands add.
    Endpoints are ε = 0 so ``SpectrumSpec`` interpolation looks like a
    scanned UV-Vis trace rather than a single spike.
    """
    if not bands:
        raise ValueError("need at least one (nm, epsilon) band")

    parsed = []
    for band in bands:
        if len(band) == 2:
            lam, eps = band
            width = fwhm_cm
        elif len(band) == 3:
            lam, eps, width = band
        else:
            raise ValueError("band must be (nm, eps) or (nm, eps, fwhm_cm)")
        parsed.append((float(lam), float(eps), float(width)))

    def epsilon_cm(wavelength_nm):
        wavenumber = 1.0e7 / wavelength_nm
        total = 0.0
        for peak_nm, peak_eps, width in parsed:
            peak_wn = 1.0e7 / peak_nm
            total += peak_eps * math.exp(-_LN2_4 * ((wavenumber - peak_wn) / width) ** 2)
        return total

    lo = min(1.0e7 / (1.0e7 / peak_nm + 2.2 * width) for peak_nm, _, width in parsed)
    hi = max(
        1.0e7 / max(1.0e7 / peak_nm - 2.2 * width, 400.0)
        for peak_nm, _, width in parsed
    )
    lo = max(200.0, lo)
    hi = min(1500.0, hi)

    grid = set()
    sample = lo
    while sample <= hi + 0.5 * step_nm:
        grid.add(round(sample, 3))
        sample += step_nm
    for peak_nm, _, _ in parsed:
        grid.add(round(peak_nm, 3))

    max_eps = max(peak_eps for _, peak_eps, _ in parsed)
    cutoff = 0.005 * max_eps
    samples = [(wavelength, epsilon_cm(wavelength)) for wavelength in sorted(grid)]
    while len(samples) > 2 and samples[0][1] < cutoff:
        samples.pop(0)
    while len(samples) > 2 and samples[-1][1] < cutoff:
        samples.pop()

    left = max(200.0, samples[0][0] - step_nm)
    right = samples[-1][0] + step_nm
    points = [(left, 0.0), *samples, (right, 0.0)]
    return vis_nm(*points, extrapolate="flat")


def compound(formula, phase, *, mp=None, bp=None, charge=0, spectrum=None):
    return Compound(
        formula,
        phase_point_list=phase_at(phase),
        mp=mp,
        bp=bp,
        charge=charge,
        spectrum=spectrum,
    )
