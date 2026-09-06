"""Build library half-reactions with E° vs SHE and every concentration left at 0."""

from .._half_reaction import HalfReaction
from ..compounds._builders import library


def hr(spec, *, E0, name="", **kwargs):
    """``HalfReaction.from_string`` with all species concentrations 0."""
    return HalfReaction.from_string(spec, E0=E0, name=name, **kwargs)
