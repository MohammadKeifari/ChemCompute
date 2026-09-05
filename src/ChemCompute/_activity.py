"""Ionic activity coefficient models for equilibrium calculations."""

from __future__ import annotations

import math
import warnings
from typing import Optional, Union

import numpy as np

R_GAS = 8.3145
VALID_ACTIVITY_MODELS = ("debye_huckel", "davies", "pitzer_lite")


def debye_huckel_A(temperature: float) -> float:
    """Debye-Hückel A parameter for log10(gamma) in aqueous solutions."""
    t_c = temperature - 273.15
    return 0.4883 + 0.0001 * t_c - 0.0000017 * t_c ** 2


def debye_huckel_B(temperature: float) -> float:
    """Debye-Hückel B parameter (1/(cm*sqrt(mol/L)))."""
    t_c = temperature - 273.15
    return 0.3241 + 0.00015 * t_c


def resolve_charge(env, compound) -> int:
    """Return ionic charge; env.charge_map overrides Compound.charge."""
    charge_map = getattr(env, "charge_map", None) or {}
    if compound.formula in charge_map:
        return int(charge_map[compound.formula])
    return int(getattr(compound, "charge", 0) or 0)


def ionic_strength(env, concentrations: np.ndarray, temperature: float) -> float:
    """
    Compute ionic strength I = 0.5 sum(c_i z_i^2) for charged aqueous species.
    """
    total = 0.0
    for j, compound in enumerate(env.compounds):
        if resolve_charge(env, compound) == 0:
            continue
        phase = compound.phase(temperature)
        if phase not in (None, "aq"):
            continue
        c = float(concentrations[j])
        if c <= 0:
            continue
        z = resolve_charge(env, compound)
        total += c * z * z
    return 0.5 * total


def _log10_gamma_debye_huckel(z: int, ionic_i: float, temperature: float, ion_size: float = 3.0) -> float:
    if z == 0 or ionic_i <= 0:
        return 0.0
    sqrt_i = math.sqrt(ionic_i)
    a_param = debye_huckel_A(temperature)
    b_param = debye_huckel_B(temperature)
    return -a_param * z * z * sqrt_i / (1.0 + b_param * ion_size * sqrt_i)


def _log10_gamma_davies(z: int, ionic_i: float, temperature: float) -> float:
    if z == 0 or ionic_i <= 0:
        return 0.0
    sqrt_i = math.sqrt(ionic_i)
    a_param = debye_huckel_A(temperature)
    return -a_param * z * z * (sqrt_i / (1.0 + sqrt_i) - 0.3 * ionic_i)


def _log10_gamma_pitzer_lite(
    z: int,
    ionic_i: float,
    temperature: float,
    pitzer_beta0: float = 0.0765,
    pitzer_pair_charge: Optional[tuple[int, int]] = None,
    resolved_charge: Optional[int] = None,
) -> float:
    """Davies base plus a single dominant beta0 term for a 1:1 electrolyte pair."""
    log_gamma = _log10_gamma_davies(z, ionic_i, temperature)
    if ionic_i <= 0 or pitzer_beta0 == 0:
        return log_gamma
    if pitzer_pair_charge is not None and resolved_charge is not None:
        z1, z2 = pitzer_pair_charge
        if abs(resolved_charge) in (abs(z1), abs(z2)):
            log_gamma -= pitzer_beta0 * ionic_i
    return log_gamma


class ActivityModel:
    """
    Activity coefficient calculator for an environment.

    Parameters
    ----------
    model : str
        ``debye_huckel``, ``davies``, or ``pitzer_lite``.
    ion_size : float
        Ion-size parameter for extended Debye-Hückel (Å, default 3.0).
    pitzer_beta0 : float
        Pitzer beta0 coefficient for pitzer_lite (default NaCl-like 0.0765).
    pitzer_pair_charge : tuple[int, int], optional
        Dominant cation/anion charges for pitzer_lite, e.g. (1, -1).
    warn_high_ionic_strength : float
        Warn when I exceeds this value (default 0.5 M).
    """

    def __init__(
        self,
        model: str = "davies",
        *,
        ion_size: float = 3.0,
        pitzer_beta0: float = 0.0765,
        pitzer_pair_charge: Optional[tuple[int, int]] = None,
        warn_high_ionic_strength: float = 0.5,
    ):
        if model not in VALID_ACTIVITY_MODELS:
            raise ValueError(f"Unknown activity model: {model}. Choose from {VALID_ACTIVITY_MODELS}.")
        self.model = model
        self.ion_size = ion_size
        self.pitzer_beta0 = pitzer_beta0
        self.pitzer_pair_charge = pitzer_pair_charge or (1, -1)
        self.warn_high_ionic_strength = warn_high_ionic_strength
        self._warned_high_i = False

    def gamma_array(self, env, concentrations: np.ndarray, temperature: float) -> np.ndarray:
        """Return activity coefficients (dimensionless) for each compound."""
        conc = np.asarray(concentrations, dtype=float)
        ionic_i = ionic_strength(env, conc, temperature)
        if (
            not self._warned_high_i
            and self.warn_high_ionic_strength is not None
            and ionic_i > self.warn_high_ionic_strength
            and self.model != "davies"
        ):
            warnings.warn(
                f"Ionic strength I={ionic_i:.3g} M exceeds {self.warn_high_ionic_strength} M; "
                "consider davies or verifying model validity.",
                stacklevel=2,
            )
            self._warned_high_i = True

        gammas = np.ones(len(env.compounds), dtype=float)
        for j, compound in enumerate(env.compounds):
            z = resolve_charge(env, compound)
            if z == 0:
                continue
            phase = compound.phase(temperature)
            if phase not in (None, "aq"):
                continue
            if self.model == "debye_huckel":
                log10_g = _log10_gamma_debye_huckel(z, ionic_i, temperature, self.ion_size)
            elif self.model == "davies":
                log10_g = _log10_gamma_davies(z, ionic_i, temperature)
            else:
                log10_g = _log10_gamma_pitzer_lite(
                    z,
                    ionic_i,
                    temperature,
                    self.pitzer_beta0,
                    self.pitzer_pair_charge,
                    z,
                )
            gammas[j] = 10.0 ** log10_g
        return gammas

    def effective_concentrations(
        self, env, concentrations: np.ndarray, temperature: float
    ) -> np.ndarray:
        """Return gamma * c for aqueous species; unchanged for others."""
        conc = np.asarray(concentrations, dtype=float)
        gammas = self.gamma_array(env, conc, temperature)
        return conc * gammas


def normalize_activity_model(
    activity_model: Optional[Union[str, ActivityModel, bool]],
) -> Optional[ActivityModel]:
    """Convert user activity_model argument to ActivityModel or None."""
    if activity_model is None or activity_model is False:
        return None
    if activity_model is True:
        return ActivityModel("davies")
    if isinstance(activity_model, ActivityModel):
        return activity_model
    if isinstance(activity_model, str):
        return ActivityModel(activity_model)
    raise TypeError("activity_model must be None, str, ActivityModel, or bool.")
