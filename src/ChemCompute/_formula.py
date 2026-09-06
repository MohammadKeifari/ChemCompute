"""Formula normalization and ionic charge inference for reaction string parsing."""

from __future__ import annotations

import re
from typing import Optional

_PHASE_SUFFIXES = (".aq", ".s", ".l", ".g")
_CHARGE_SUFFIX_POS = re.compile(r"\+(\d*)$")
_CHARGE_SUFFIX_NEG = re.compile(r"-(\d*)$")


def split_phase_suffix(name: str) -> tuple[str, Optional[str]]:
    """Return ``(core, phase)`` where phase is one of s/l/g/aq or None."""
    for suffix, phase in ((".aq", "aq"), (".s", "s"), (".l", "l"), (".g", "g")):
        if name.endswith(suffix):
            return name[: -len(suffix)], phase
    return name, None


def normalize_formula(name: str) -> str:
    """
    Canonicalize a species token: strip outer brackets, keep charge suffix in name.

    ``[Fe(CN)6]-4`` and ``Fe(CN)6-4`` both become ``Fe(CN)6-4``.
    """
    core, phase = split_phase_suffix(name.strip())
    bracket_with_charge = re.match(r"^\[(.+)\]([+-]\d*)$", core)
    if bracket_with_charge:
        core = bracket_with_charge.group(1) + bracket_with_charge.group(2)
    elif core.startswith("[") and core.endswith("]"):
        core = core[1:-1]
    if phase is not None:
        core = f"{core}.{phase}"
    return core


def infer_ionic_charge(formula: str) -> int:
    """
    Infer signed ionic charge from trailing ``+`` / ``-`` notation.

    Examples: ``H+`` → 1, ``OH-`` → -1, ``Fe+3`` → 3, ``SeO4-2`` → -2,
    ``Fe(CN)6-4`` → -4. Neutral species return 0.
    """
    core, _ = split_phase_suffix(formula.strip())
    bracket_with_charge = re.match(r"^\[(.+)\]([+-]\d*)$", core)
    if bracket_with_charge:
        core = bracket_with_charge.group(1) + bracket_with_charge.group(2)
    elif core.startswith("[") and core.endswith("]"):
        core = core[1:-1]

    match = _CHARGE_SUFFIX_POS.search(core)
    if match:
        digits = match.group(1)
        return int(digits) if digits else 1

    match = _CHARGE_SUFFIX_NEG.search(core)
    if match:
        digits = match.group(1)
        return -int(digits) if digits else -1

    return 0


def compound_from_species_token(name: str, *, T: float = 298):
    """Build a :class:`Compound` with inferred charge and phase from a parsed token."""
    from ._general import Compound

    canonical = normalize_formula(name)
    core, phase = split_phase_suffix(canonical)
    charge = infer_ionic_charge(core)

    phase_point_list = None
    if phase is not None:
        phase_point_list = [{"temperature": T, "phase": phase}]

    return Compound(formula=core, phase_point_list=phase_point_list, charge=charge)
