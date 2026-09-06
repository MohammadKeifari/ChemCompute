from ._general import Compound, Enviroment, Reaction, XS
from ._equilibrium import EquilibriumResult
from ._half_reaction import HalfReaction, BoundaryLine
from ._pourbaix import (
    Pourbaix,
    PourbaixFrameIntersection,
    PourbaixMetadata,
    PourbaixResult,
    infer_pourbaix_metadata,
)
from ._pourbaix_graph import (
    PourbaixBoundary,
    PourbaixGraph,
    PourbaixJunction,
    boundary_Eh,
    build_pourbaix_graph,
    format_junction_label,
    format_junction_plot_label,
    graph_speciation,
)
from ._activity import ActivityModel, VALID_ACTIVITY_MODELS, ionic_strength, resolve_charge
from ._buffer import BufferDiagnostics, BufferPairDiagnostics, buffer_diagnostics
from ._mixing import ScaledEnviroment, combine_environments
from ._titration import Titration, TitrationResult, mix_sample_with_titrant
from ._uvvis import SpectrumSpec, uvvis_spectrum
from ._bio_kinetics import integrate_bio_kinetics, uses_bio_kinetics, RATE_LAW_FUNCTIONS
from . import bio_templates
from . import compounds
from . import environments
from . import half_reactions
from . import reactions
from .bio_templates import (
    single_substrate_mm,
    competitive_inhibition,
    uncompetitive_inhibition,
    noncompetitive_inhibition,
    mixed_inhibition,
    sequential_pathway,
)

Environment = Enviroment

__all__ = [
    "ActivityModel",
    "BoundaryLine",
    "BufferDiagnostics",
    "BufferPairDiagnostics",
    "Compound",
    "Enviroment",
    "Environment",
    "EquilibriumResult",
    "HalfReaction",
    "Pourbaix",
    "PourbaixBoundary",
    "PourbaixFrameIntersection",
    "PourbaixGraph",
    "PourbaixJunction",
    "PourbaixMetadata",
    "PourbaixResult",
    "RATE_LAW_FUNCTIONS",
    "Reaction",
    "ScaledEnviroment",
    "SpectrumSpec",
    "Titration",
    "TitrationResult",
    "VALID_ACTIVITY_MODELS",
    "XS",
    "bio_templates",
    "boundary_Eh",
    "buffer_diagnostics",
    "build_pourbaix_graph",
    "combine_environments",
    "compounds",
    "competitive_inhibition",
    "environments",
    "format_junction_label",
    "format_junction_plot_label",
    "graph_speciation",
    "half_reactions",
    "infer_pourbaix_metadata",
    "integrate_bio_kinetics",
    "ionic_strength",
    "mix_sample_with_titrant",
    "mixed_inhibition",
    "noncompetitive_inhibition",
    "reactions",
    "resolve_charge",
    "sequential_pathway",
    "single_substrate_mm",
    "uncompetitive_inhibition",
    "uses_bio_kinetics",
    "uvvis_spectrum",
]
