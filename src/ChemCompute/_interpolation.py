"""Live Compound slots for Reaction.from_string / HalfReaction.from_string."""

from __future__ import annotations

import re
import weakref

_COMPOUNDS_BY_ID: weakref.WeakValueDictionary[int, object] = weakref.WeakValueDictionary()
_SLOT = re.compile(r"^@c(\d+)$")
_SLOT_WITH_PHASE = re.compile(r"^@c\d+\.(s|l|g|aq)$")
_LIVE_SECTION = re.compile(
    r"^(\d+(?:\.\d+)?_)?@c\d+(_-?\d+(?:\.\d+)?)?(\.s|\.g|\.l|\.aq)?$"
)


def register_compound(compound) -> str:
    """Store ``compound`` and return an ``@c{id}`` slot for reaction strings."""
    _COMPOUNDS_BY_ID[id(compound)] = compound
    return f"@c{id(compound)}"


def resolve_compound_token(name: str):
    """Return the live Compound for ``@c123``, or None if ``name`` is a normal formula."""
    match = _SLOT.match(name)
    if match is None:
        return None
    compound = _COMPOUNDS_BY_ID.get(int(match.group(1)))
    if compound is None:
        raise ValueError(f"Unknown compound token {name!r}")
    return compound


def is_live_compound_section(section: str) -> bool:
    """True when a parsed species term is (or contains) an ``@cN`` slot."""
    return bool(_LIVE_SECTION.match(section))


def compound_from_parsed_name(name: str, *, T: float = 298):
    """Resolve an ``@cN`` slot or build a Compound from a species token."""
    from ._formula import compound_from_species_token
    from ._half_reaction import _reject_electron_formula

    if _SLOT_WITH_PHASE.match(name):
        raise ValueError(
            "Phase suffix is not allowed on interpolated compounds; "
            "phase lives on the Compound object."
        )
    live = resolve_compound_token(name)
    if live is not None:
        return live
    _reject_electron_formula(name)
    return compound_from_species_token(name, T=T)
