"""Aqueous complex redox couples vs SHE (25 °C)."""

from ..compounds import fe_cn6_3minus, fe_cn6_4minus
from ._builders import hr, library


# Fe(CN)6-3 + e- ⇌ Fe(CN)6-4.
@library
def ferricyanide():
    return hr(
        f"{fe_cn6_3minus().token} & @e = {fe_cn6_4minus().token}",
        E0=0.356,
        name="Fe(CN)6-3/Fe(CN)6-4",
    )


__all__ = [
    "ferricyanide",
]
