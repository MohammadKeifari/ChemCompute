"""Build an Enviroment that uses library reactions and compounds as-is."""

from __future__ import annotations

from .._general import Enviroment


def assemble(*library_items, concentrations=None, T=298, volume=1.0):
    """Return an environment wired to the given library reactions or half-reactions (no copies)."""
    resolved = []
    for item in library_items:
        if callable(item):
            item = item()
        resolved.append(item)
    kwargs = {"T": T, "volume": volume}
    if concentrations:
        kwargs["concentrations"] = concentrations
    return Enviroment(*resolved, **kwargs)
