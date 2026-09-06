"""Coordination complexes and acid–base indicators with UV-Vis envelopes."""

from ._builders import compound, vis_bands, vis_nm

# Fe(III) thiocyanate — Frank & Oswalt 447 nm; broad LMCT toward ~470 nm
fescn = compound(
    "FeSCN+2", "aq", charge=2, spectrum=vis_bands((447, 4700.0, 4500.0))
)

# hexacyanoferrates
fe_cn6_3minus = compound(
    "Fe(CN)6-3",
    "aq",
    charge=-3,
    spectrum=vis_bands((260, 1200.0, 4000.0), (420, 1040.0, 3000.0)),
)
fe_cn6_4minus = compound("Fe(CN)6-4", "aq", charge=-4)

# ammine complexes
cu_nh3_4_2plus = compound(
    "Cu(NH3)4+2", "aq", charge=2, spectrum=vis_bands((600, 56.0), fwhm_cm=3500.0)
)
ni_nh3_6_2plus = compound(
    "Ni(NH3)6+2", "aq", charge=2, spectrum=vis_bands((355, 6.5), (571, 4.8))
)
co_nh3_6_3plus = compound(
    "Co(NH3)6+3", "aq", charge=3, spectrum=vis_bands((339, 46.0), (475, 56.0))
)
cr_nh3_6_3plus = compound(
    "Cr(NH3)6+3", "aq", charge=3, spectrum=vis_bands((350, 33.0), (462, 40.0))
)
ag_nh3_2_plus = compound("Ag(NH3)2+", "aq", charge=1)
zn_nh3_4_2plus = compound("Zn(NH3)4+2", "aq", charge=2)

# tetrahalometallates
co_cl4_2minus = compound(
    "CoCl4-2",
    "aq",
    charge=-2,
    spectrum=vis_nm(
        (550, 0.0),
        (590, 150.0),
        (625, 380.0),
        (655, 500.0),
        (670, 560.0),
        (690, 600.0),
        (720, 180.0),
        (760, 40.0),
        (800, 0.0),
    ),
)
ni_cl4_2minus = compound("NiCl4-2", "aq", charge=-2)
cu_cl4_2minus = compound(
    "CuCl4-2", "aq", charge=-2, spectrum=vis_bands((400, 80.0, 4000.0))
)

# phenanthroline / bipyridine
ferroin = compound(
    "Fe(phen)3+2",
    "aq",
    charge=2,
    spectrum=vis_nm(
        (380, 0.0),
        (420, 2500.0),
        (450, 5500.0),
        (475, 8500.0),
        (510, 11100.0),
        (540, 5000.0),
        (580, 800.0),
        (640, 0.0),
    ),
)
ru_bpy3_2plus = compound(
    "Ru(bpy)3+2",
    "aq",
    charge=2,
    spectrum=vis_bands((285, 87000.0, 3000.0), (452, 14600.0, 2000.0)),
)

# oxalate / EDTA-type teaching complexes
fe_ox3_3minus = compound("Fe(C2O4)3-3", "aq", charge=-3)
ni_en3_2plus = compound(
    "Ni(en)3+2", "aq", charge=2, spectrum=vis_bands((545, 6.5))
)
ferrocene = compound(
    "Fe(C5H5)2",
    "s",
    mp=446.0,
    bp=522.0,
    spectrum=vis_bands((325, 51.0, 3000.0), (440, 91.0, 2800.0)),
)

# indicators / dyes (aqueous envelopes around published λmax)
methylene_blue = compound(
    "MB+",
    "aq",
    charge=1,
    spectrum=vis_nm(
        (520, 0.0),
        (560, 8000.0),
        (590, 22000.0),
        (610, 38000.0),
        (640, 48000.0),
        (664, 74000.0),
        (685, 22000.0),
        (710, 3000.0),
        (750, 0.0),
    ),
)
crystal_violet = compound(
    "CV+", "aq", charge=1, spectrum=vis_bands((590, 87000.0, 1200.0))
)
fluorescein = compound(
    "Fl-2", "aq", charge=-2, spectrum=vis_bands((490, 76900.0, 1400.0))
)
phenolphthalein_pink = compound(
    "HIn-2", "aq", charge=-2, spectrum=vis_bands((552, 26000.0, 1500.0))
)
methyl_orange = compound(
    "MO-", "aq", charge=-1, spectrum=vis_bands((464, 22600.0, 1800.0))
)
bromothymol_blue = compound(
    "BTB-", "aq", charge=-1, spectrum=vis_bands((616, 32500.0, 1400.0))
)

__all__ = [
    "ag_nh3_2_plus",
    "bromothymol_blue",
    "co_cl4_2minus",
    "co_nh3_6_3plus",
    "cr_nh3_6_3plus",
    "crystal_violet",
    "cu_cl4_2minus",
    "cu_nh3_4_2plus",
    "fe_cn6_3minus",
    "fe_cn6_4minus",
    "fe_ox3_3minus",
    "ferroin",
    "ferrocene",
    "fescn",
    "fluorescein",
    "methyl_orange",
    "methylene_blue",
    "ni_cl4_2minus",
    "ni_en3_2plus",
    "ni_nh3_6_2plus",
    "phenolphthalein_pink",
    "ru_bpy3_2plus",
    "zn_nh3_4_2plus",
]
