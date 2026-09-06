"""Build library reactions with K set and every concentration left at 0."""

from .._general import Reaction
from ..compounds._builders import library


def rxn(spec, *, K=1.0, infinite_K=False, **kwargs):
    """``Reaction.from_string`` with an empty concentration map (all species 0)."""
    return Reaction.from_string(
        spec,
        concentrations={},
        K=K,
        infinite_K=infinite_K,
        **kwargs,
    )
