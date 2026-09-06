"""Pourbaix diagram scanning via fixed pH and electrode potential."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional, Sequence

import numpy as np

from ._general import POURBAIX_RESERVED_SPECIES
from ._half_reaction import BoundaryLine
from ._pourbaix_graph import (
    NERNST_K,
    JunctionLabelStyle,
    PourbaixBoundary,
    PourbaixGraph,
    PourbaixJunction,
    assign_junction_point_ids,
    boundary_species_at_pH,
    build_pourbaix_graph,
    compute_analytic_geometry,
    element_totals_from_env,
    format_junction_plot_label,
    graph_speciation,
)
from ._titration import _find_h_plus_index

BoundaryMode = Literal["dominant", "all"]


def _pH_grid(pH_min: float, pH_max: float, steps: int) -> np.ndarray:
    if steps < 2:
        raise ValueError("pH grid steps must be at least 2.")
    if pH_min > pH_max:
        raise ValueError("Require pH_min <= pH_max.")
    return np.linspace(pH_min, pH_max, int(steps))


def _eh_grid(eh_min: float, eh_max: float, steps: int) -> np.ndarray:
    if steps < 2:
        raise ValueError("Eh grid steps must be at least 2.")
    if eh_min > eh_max:
        raise ValueError("Require Eh_min <= Eh_max.")
    return np.linspace(eh_min, eh_max, int(steps))


def _dominant_index(concentrations: Sequence[float], env, track_species: Sequence[str]) -> int:
    label_map = {compound.formula: j for j, compound in enumerate(env.compounds)}
    best_index = 0
    best_value = -1.0
    for index, label in enumerate(track_species):
        species_index = label_map.get(label)
        if species_index is None:
            continue
        value = float(concentrations[species_index])
        if value > best_value:
            best_value = value
            best_index = index
    return best_index


def apply_pourbaix_grid_state(
    env,
    pH: float,
    Eh: float,
    overrides: Optional[dict[str, float]] = None,
) -> None:
    """Fix pH, Eh, and optional species overrides. Always owns H+ and OH-."""
    h_conc = 10.0 ** (-float(pH))
    oh_conc = 1e-14 / max(h_conc, 1e-300)
    env.set_buffer(["H+"])
    env.set_electrode_potential(Eh=float(Eh))
    if overrides:
        safe = {
            key: value
            for key, value in overrides.items()
            if (key if isinstance(key, str) else getattr(key, "formula", "")) not in POURBAIX_RESERVED_SPECIES
        }
        if safe:
            env._apply_concentration_overrides(safe)
    reserved = {"H+": h_conc}
    if "OH-" in env.compound_labels:
        reserved["OH-"] = oh_conc
    env._apply_concentration_overrides(reserved, allow_pourbaix_reserved=True)


def _neighbor_concentrations(
    speciation: np.ndarray,
    track_species: Sequence[str],
    i: int,
    j: int,
) -> Optional[dict[str, float]]:
    for ni, nj in ((i, j - 1), (i - 1, j), (i - 1, j - 1)):
        if ni < 0 or nj < 0:
            continue
        if np.all(speciation[ni, nj] <= 0):
            continue
        return {
            label: float(speciation[ni, nj, index])
            for index, label in enumerate(track_species)
        }
    return None


def _adjacent_species_pairs(grid_dominant: np.ndarray) -> set[tuple[int, int]]:
    pairs: set[tuple[int, int]] = set()
    n_pH, n_eh = grid_dominant.shape
    for i in range(n_pH):
        for j in range(n_eh):
            current = int(grid_dominant[i, j])
            if i + 1 < n_pH:
                other = int(grid_dominant[i + 1, j])
                if other != current:
                    pairs.add((min(current, other), max(current, other)))
            if j + 1 < n_eh:
                other = int(grid_dominant[i, j + 1])
                if other != current:
                    pairs.add((min(current, other), max(current, other)))
    return pairs


def _split_active_segments(
    pH: np.ndarray,
    Eh: np.ndarray,
    active: np.ndarray,
) -> list[tuple[np.ndarray, np.ndarray]]:
    segments: list[tuple[np.ndarray, np.ndarray]] = []
    start: Optional[int] = None
    for index, include in enumerate(active):
        if include and start is None:
            start = index
        elif not include and start is not None:
            segments.append((pH[start:index], Eh[start:index]))
            start = None
    if start is not None:
        segments.append((pH[start:], Eh[start:]))
    return segments


def _on_dominant_interface(
    pH: float,
    Eh: float,
    idx_a: int,
    idx_b: int,
    graph: PourbaixGraph,
    totals: dict[str, float],
    track_species: Sequence[str],
    *,
    background_ions: Optional[dict[str, float]] = None,
) -> bool:
    """True when (pH, Eh) lies on the border between two dominant region species."""
    label_to_index = {name: index for index, name in enumerate(track_species)}

    def dominant_index(ph: float, eh: float) -> int:
        _, index = graph_speciation(
            ph,
            eh,
            graph,
            totals,
            background_ions=background_ions,
        )
        name = graph.all_track_species[index]
        return label_to_index.get(name, -1)

    seen = {
        dominant_index(pH, Eh),
        dominant_index(pH, Eh + 1e-3),
        dominant_index(pH, Eh - 1e-3),
        dominant_index(pH + 1e-3, Eh),
        dominant_index(pH - 1e-3, Eh),
    }
    seen.discard(-1)
    return idx_a in seen and idx_b in seen


def _clip_curve_to_dominant_pair(
    pH_curve: np.ndarray,
    eh_curve: np.ndarray,
    idx_a: int,
    idx_b: int,
    name_a: str,
    name_b: str,
    boundary: PourbaixBoundary,
    graph: PourbaixGraph,
    totals: dict[str, float],
    track_species: Sequence[str],
    *,
    pH_min: float,
    pH_max: float,
    eh_min: float,
    eh_max: float,
    background_ions: Optional[dict[str, float]] = None,
) -> list[tuple[np.ndarray, np.ndarray]]:
    active = np.zeros(len(pH_curve), dtype=bool)
    target = {name_a, name_b}
    for index in range(len(pH_curve)):
        pH = float(pH_curve[index])
        Eh = float(eh_curve[index])
        if pH < pH_min or pH > pH_max or Eh < eh_min or Eh > eh_max:
            continue
        left, right = boundary_species_at_pH(boundary, pH, graph)
        if {left, right} != target:
            continue
        if not _on_dominant_interface(
            pH,
            Eh,
            idx_a,
            idx_b,
            graph,
            totals,
            track_species,
            background_ions=background_ions,
        ):
            continue
        active[index] = True
    return _split_active_segments(pH_curve, eh_curve, active)


def _dominant_analytic_boundary_lines(
    grid_dominant: np.ndarray,
    track_species: Sequence[str],
    analytic_boundaries: Sequence[PourbaixBoundary],
    graph: PourbaixGraph,
    totals: dict[str, float],
    *,
    pH_min: float,
    pH_max: float,
    eh_min: float,
    eh_max: float,
    background_ions: Optional[dict[str, float]] = None,
) -> list[tuple[int, int, np.ndarray, np.ndarray]]:
    """
    Smooth analytic equilibrium curves, clipped to segments where the pair
    separates neighboring predominant regions on the diagram.
    """
    lines: list[tuple[int, int, np.ndarray, np.ndarray]] = []
    for idx_a, idx_b in sorted(_adjacent_species_pairs(grid_dominant)):
        name_a = track_species[idx_a]
        name_b = track_species[idx_b]
        for boundary in analytic_boundaries:
            if boundary.kind == "water":
                continue
            for pH_seg, eh_seg in _clip_curve_to_dominant_pair(
                boundary.pH,
                boundary.Eh,
                idx_a,
                idx_b,
                name_a,
                name_b,
                boundary,
                graph,
                totals,
                track_species,
                pH_min=pH_min,
                pH_max=pH_max,
                eh_min=eh_min,
                eh_max=eh_max,
                background_ions=background_ions,
            ):
                if len(pH_seg) >= 2:
                    lines.append((idx_a, idx_b, pH_seg, eh_seg))
    return lines


def _contour_boundary_lines_for_pair(
    grid_pH: np.ndarray,
    grid_Eh: np.ndarray,
    speciation: np.ndarray,
    grid_dominant: np.ndarray,
    track_species: Sequence[str],
    left: int,
    right: int,
    graph: PourbaixGraph,
    totals: dict[str, float],
    *,
    pH_min: float,
    pH_max: float,
    eh_min: float,
    eh_max: float,
    background_ions: Optional[dict[str, float]] = None,
) -> list[tuple[int, int, np.ndarray, np.ndarray]]:
    pH_mesh, eh_mesh = np.meshgrid(grid_pH, grid_Eh, indexing="ij")
    diff = speciation[:, :, left] - speciation[:, :, right]
    lines: list[tuple[int, int, np.ndarray, np.ndarray]] = []
    try:
        import matplotlib.pyplot as plt

        fig = plt.figure()
        try:
            cs = plt.contour(pH_mesh, eh_mesh, diff, levels=[0.0])
            for collection in cs.collections:
                for path in collection.get_paths():
                    vertices = path.vertices
                    if len(vertices) < 2:
                        continue
                    active = np.array(
                        [
                            _on_dominant_interface(
                                float(vertices[k, 0]),
                                float(vertices[k, 1]),
                                left,
                                right,
                                graph,
                                totals,
                                track_species,
                                background_ions=background_ions,
                            )
                            for k in range(len(vertices))
                        ],
                        dtype=bool,
                    )
                    for pH_seg, eh_seg in _split_active_segments(vertices[:, 0], vertices[:, 1], active):
                        if len(pH_seg) >= 2:
                            lines.append((left, right, pH_seg, eh_seg))
        finally:
            plt.close(fig)
    except Exception:
        return []
    return lines


def _equal_concentration_boundaries(
    grid_pH: np.ndarray,
    grid_Eh: np.ndarray,
    speciation: np.ndarray,
    grid_dominant: np.ndarray,
    track_species: Sequence[str],
    *,
    analytic_boundaries: Optional[Sequence[PourbaixBoundary]] = None,
    graph: Optional[PourbaixGraph] = None,
    totals: Optional[dict[str, float]] = None,
    pH_min: float = 0.0,
    pH_max: float = 14.0,
    eh_min: float = -1.5,
    eh_max: float = 1.5,
    background_ions: Optional[dict[str, float]] = None,
) -> list[tuple[int, int, np.ndarray, np.ndarray]]:
    """Analytic point-to-point lines between dominant neighbors only."""
    if analytic_boundaries and graph is not None and totals is not None:
        lines = _dominant_analytic_boundary_lines(
            grid_dominant,
            track_species,
            analytic_boundaries,
            graph,
            totals,
            pH_min=pH_min,
            pH_max=pH_max,
            eh_min=eh_min,
            eh_max=eh_max,
            background_ions=background_ions,
        )
        if lines:
            return lines

    if graph is None or totals is None:
        return []

    lines: list[tuple[int, int, np.ndarray, np.ndarray]] = []
    for left, right in sorted(_adjacent_species_pairs(grid_dominant)):
        lines.extend(
            _contour_boundary_lines_for_pair(
                grid_pH,
                grid_Eh,
                speciation,
                grid_dominant,
                track_species,
                left,
                right,
                graph,
                totals,
                pH_min=pH_min,
                pH_max=pH_max,
                eh_min=eh_min,
                eh_max=eh_max,
                background_ions=background_ions,
            )
        )
    return lines


def _dominant_junction_points(
    grid_pH: np.ndarray,
    grid_Eh: np.ndarray,
    grid_dominant: np.ndarray,
    track_species: Sequence[str],
    *,
    junction_labels: Optional[dict[tuple[str, ...], str]] = None,
) -> list[PourbaixJunction]:
    """
    Triple points where three or more predominant regions meet on the diagram grid.
    """
    from ._pourbaix_graph import assign_junction_point_ids, format_junction_label

    raw: list[tuple[float, float, tuple[str, ...]]] = []
    n_pH, n_eh = grid_dominant.shape
    for i in range(n_pH - 1):
        for j in range(n_eh - 1):
            indices = {
                int(grid_dominant[i, j]),
                int(grid_dominant[i + 1, j]),
                int(grid_dominant[i, j + 1]),
                int(grid_dominant[i + 1, j + 1]),
            }
            if len(indices) < 3:
                continue
            pH = 0.5 * (float(grid_pH[i]) + float(grid_pH[i + 1]))
            Eh = 0.5 * (float(grid_Eh[j]) + float(grid_Eh[j + 1]))
            species = tuple(sorted(track_species[index] for index in indices))
            raw.append((pH, Eh, species))

    if not raw:
        return []

    clusters: list[tuple[float, float, set[str]]] = []
    eps_pH = max(float(grid_pH[1] - grid_pH[0]), 1e-6) * 1.5
    eps_eh = max(float(grid_Eh[1] - grid_Eh[0]), 1e-6) * 1.5
    for pH, Eh, species_tuple in raw:
        merged = False
        for index, (cx, cy, species_set) in enumerate(clusters):
            if abs(pH - cx) <= eps_pH and abs(Eh - cy) <= eps_eh:
                clusters[index] = (
                    (cx + pH) / 2.0,
                    (cy + Eh) / 2.0,
                    species_set | set(species_tuple),
                )
                merged = True
                break
        if not merged:
            clusters.append((pH, Eh, set(species_tuple)))

    junctions: list[PourbaixJunction] = []
    for pH, Eh, species_set in clusters:
        if len(species_set) < 3:
            continue
        species_tuple = tuple(sorted(species_set))
        label = format_junction_label(species_tuple)
        if junction_labels and species_tuple in junction_labels:
            label = junction_labels[species_tuple]
        junctions.append(
            PourbaixJunction(
                pH=pH,
                Eh=Eh,
                species=species_tuple,
                label=label,
            )
        )
    return assign_junction_point_ids(junctions)


@dataclass
class PourbaixMetadata:
    """Legacy metadata view derived from :class:`PourbaixGraph`."""

    track_species: list[str]
    oxidation_groups: list[list[str]]
    redox_chains: list[list[str]]
    pka_pairs: list[tuple[str, str, float]]
    graph: PourbaixGraph = field(repr=False)


def infer_pourbaix_metadata(env) -> PourbaixMetadata:
    """
    Derive Pourbaix bookkeeping from an environment's half-reactions and
    acid-base reactions. Prefer :func:`build_pourbaix_graph` for new code.
    """
    graph = build_pourbaix_graph(env)
    if not graph.chains:
        raise ValueError("Pourbaix inference requires at least one half-reaction.")
    oxidation_groups: list[list[str]] = []
    redox_chains: list[list[str]] = []
    for chain in graph.chains:
        oxidation_groups.extend(chain.oxidation_levels)
        redox_chains.extend(chain.oxidation_levels)
    return PourbaixMetadata(
        track_species=list(graph.all_track_species),
        oxidation_groups=oxidation_groups,
        redox_chains=redox_chains,
        pka_pairs=[(e.acid, e.base, e.pka) for e in graph.acid_base_edges],
        graph=graph,
    )


@dataclass
class PourbaixResult:
    grid_pH: np.ndarray
    grid_Eh: np.ndarray
    grid_dominant: np.ndarray
    track_species: list[str]
    speciation: Optional[np.ndarray] = None
    converged: Optional[np.ndarray] = None
    geometry_source: Literal["analytic", "grid"] = "analytic"
    analytic_boundaries: list[PourbaixBoundary] = field(default_factory=list)
    junction_points: list[PourbaixJunction] = field(default_factory=list)
    analytic_junction_points: list[PourbaixJunction] = field(default_factory=list)
    equal_boundary_lines: list[tuple[int, int, np.ndarray, np.ndarray]] = field(default_factory=list)
    graph: Optional[PourbaixGraph] = None
    boundary_lines: list[BoundaryLine] = field(default_factory=list)

    def matrix(self) -> np.ndarray:
        return self.grid_dominant

    def region_matrix(self) -> np.ndarray:
        """Species index per (pH, Eh) grid cell."""
        return self.grid_dominant

    def dominant_species_at(self, pH: float, Eh: float) -> str:
        i = int(np.argmin(np.abs(self.grid_pH - pH)))
        j = int(np.argmin(np.abs(self.grid_Eh - Eh)))
        return self.track_species[int(self.grid_dominant[i, j])]

    def junction_table(self, *, source: Literal["dominant", "analytic"] = "dominant") -> list[dict]:
        """Junction coordinates keyed by ``P1``, ``P2``, … and meeting species."""
        points = self.junction_points if source == "dominant" else self.analytic_junction_points
        return [
            {
                "point_id": junction.point_id,
                "pH": junction.pH,
                "Eh": junction.Eh,
                "species": list(junction.species),
                "label": junction.label,
                "boundary_kinds": list(junction.boundary_kinds),
            }
            for junction in points
        ]

    def junction_coords(self, point_id: str, *, source: Literal["dominant", "analytic"] = "dominant") -> tuple[float, float]:
        """Return ``(pH, Eh)`` for a numbered junction such as ``P4``."""
        points = self.junction_points if source == "dominant" else self.analytic_junction_points
        for junction in points:
            if junction.point_id == point_id:
                return junction.pH, junction.Eh
        raise KeyError(f"Unknown junction point {point_id!r}.")

    def boundary_table(self) -> list[dict]:
        rows: list[dict] = []
        for boundary in self.analytic_boundaries:
            for pH, Eh in zip(boundary.pH, boundary.Eh):
                rows.append(
                    {
                        "left": boundary.left,
                        "right": boundary.right,
                        "kind": boundary.kind,
                        "pH": float(pH),
                        "Eh": float(Eh),
                    }
                )
        for junction in self.junction_points:
            rows.append(
                {
                    "kind": "junction",
                    "point_id": junction.point_id,
                    "label": junction.label,
                    "species": list(junction.species),
                    "pH": junction.pH,
                    "Eh": junction.Eh,
                }
            )
        return rows

    def _junction_annotation(
        self,
        junction: PourbaixJunction,
        *,
        style: JunctionLabelStyle,
        pH_decimals: int,
        eh_decimals: int,
    ) -> str:
        return format_junction_plot_label(
            junction,
            style=style,
            pH_decimals=pH_decimals,
            eh_decimals=eh_decimals,
        )

    def plot(
        self,
        *,
        ax=None,
        show: bool = True,
        save: Optional[str] = None,
        include_water_lines: bool = True,
        boundary_mode: BoundaryMode = "dominant",
        show_analytic_boundaries: Optional[bool] = None,
        show_equal_boundaries: Optional[bool] = None,
        show_junction_labels: bool = True,
        junction_label_style: JunctionLabelStyle = "numbered_coords",
        junction_pH_decimals: int = 2,
        junction_eh_decimals: int = 2,
        labels: Optional[dict[int, str]] = None,
    ):
        import matplotlib.pyplot as plt
        from matplotlib.colors import BoundaryNorm, ListedColormap

        if show_analytic_boundaries is None:
            show_analytic_boundaries = boundary_mode == "all"
        if show_equal_boundaries is None:
            show_equal_boundaries = boundary_mode == "dominant"

        junctions = (
            self.analytic_junction_points
            if boundary_mode == "all"
            else self.junction_points
        )

        if ax is None:
            _, ax = plt.subplots(figsize=(8, 7))

        pH_mesh, eh_mesh = np.meshgrid(self.grid_pH, self.grid_Eh, indexing="ij")
        label_map = labels or {i: name for i, name in enumerate(self.track_species)}
        n_species = max(len(self.track_species), int(self.grid_dominant.max()) + 1)
        cmap = ListedColormap(plt.cm.tab10.colors[: max(n_species, 1)])
        norm = BoundaryNorm(np.arange(-0.5, n_species + 0.5, 1), n_species)

        ax.pcolormesh(pH_mesh, eh_mesh, self.grid_dominant, cmap=cmap, norm=norm, shading="auto")

        if show_analytic_boundaries and self.analytic_boundaries:
            for boundary in self.analytic_boundaries:
                if boundary.kind == "water" and not include_water_lines:
                    continue
                ax.plot(boundary.pH, boundary.Eh, "k-", lw=0.8, alpha=0.55)

        if show_equal_boundaries and self.equal_boundary_lines:
            for _, _, pH_line, eh_line in self.equal_boundary_lines:
                ax.plot(pH_line, eh_line, "k-", lw=0.9, alpha=0.85)

        if show_junction_labels and junctions:
            for junction in junctions:
                ax.plot(junction.pH, junction.Eh, "ko", ms=4)
                ax.annotate(
                    self._junction_annotation(
                        junction,
                        style=junction_label_style,
                        pH_decimals=junction_pH_decimals,
                        eh_decimals=junction_eh_decimals,
                    ),
                    (junction.pH, junction.Eh),
                    fontsize=6,
                    xytext=(3, 3),
                    textcoords="offset points",
                )

        if include_water_lines and not (
            show_analytic_boundaries and any(b.kind == "water" for b in self.analytic_boundaries)
        ):
            ax.plot(self.grid_pH, 1.229 - 2 * NERNST_K * self.grid_pH, "k--", lw=0.8)
            ax.plot(self.grid_pH, -2 * NERNST_K * self.grid_pH, "k--", lw=0.8)

        ax.set_xlim(self.grid_pH[0], self.grid_pH[-1])
        ax.set_ylim(self.grid_Eh[0], self.grid_Eh[-1])
        ax.set_xlabel("pH")
        ax.set_ylabel("Eh (V vs SHE)")
        ax.set_title("Pourbaix diagram")
        handles = [
            plt.Line2D([0], [0], marker="s", ls="", color=cmap(i))
            for i in range(len(self.track_species))
        ]
        ax.legend(handles, [label_map[i] for i in range(len(self.track_species))], loc="upper right", fontsize=7, ncol=2)
        if save:
            plt.savefig(save, bbox_inches="tight", dpi=150)
        if show and not save:
            plt.show()
        return ax

    def plot_predominance(self, *, ax=None, show: bool = True, save: Optional[str] = None, labels=None):
        return self.plot(
            ax=ax,
            show=show,
            save=save,
            labels=labels,
            show_analytic_boundaries=False,
            show_equal_boundaries=False,
            show_junction_labels=False,
        )

    def plot_boundaries(
        self,
        *,
        ax=None,
        show: bool = True,
        save: Optional[str] = None,
        include_water_lines: bool = True,
        boundary_mode: BoundaryMode = "dominant",
        junction_label_style: JunctionLabelStyle = "numbered_coords",
        junction_pH_decimals: int = 2,
        junction_eh_decimals: int = 2,
    ):
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(8, 6))
        if boundary_mode == "all":
            for boundary in self.analytic_boundaries:
                if boundary.kind == "water" and not include_water_lines:
                    continue
                ax.plot(boundary.pH, boundary.Eh, label=f"{boundary.left} / {boundary.right}")
            junctions = self.analytic_junction_points
        else:
            for left, right, pH_line, eh_line in self.equal_boundary_lines:
                ax.plot(
                    pH_line,
                    eh_line,
                    label=f"{self.track_species[left]} / {self.track_species[right]}",
                )
            junctions = self.junction_points
        for junction in junctions:
            ax.plot(junction.pH, junction.Eh, "ko", ms=4)
            ax.annotate(
                self._junction_annotation(
                    junction,
                    style=junction_label_style,
                    pH_decimals=junction_pH_decimals,
                    eh_decimals=junction_eh_decimals,
                ),
                (junction.pH, junction.Eh),
                fontsize=7,
            )
        if include_water_lines:
            ax.plot(self.grid_pH, -NERNST_K * self.grid_pH, "--", color="gray", label="H+/H2")
            ax.plot(self.grid_pH, 1.229 - NERNST_K * self.grid_pH, "--", color="black", label="O2/H2O")
        ax.set_xlabel("pH")
        ax.set_ylabel("Eh (V vs SHE)")
        ax.set_title("Pourbaix boundaries")
        ax.legend(loc="best", fontsize=8)
        if save:
            plt.savefig(save, bbox_inches="tight")
        if show and not save:
            plt.show()
        return ax


@dataclass
class Pourbaix:
    """
    Pourbaix diagram scanner over pH and Eh.

    Builds a :class:`PourbaixGraph` from the environment's half-reactions and
    classified ``Reaction`` edges (pKa, Ksp, association, hydration).

    ``H+`` and ``OH-`` are reserved: Pourbaix sets them from pH at each grid point.

    Grid resolution
    ---------------
    ``pH_steps`` and ``Eh_steps`` set the **predominance-region grid** used for the
    colored background in **both** ``model`` and ``equilibrium`` modes. Each cell
    stores one dominant species index. Boundaries in ``model`` mode are computed
    analytically (``geometry_pH_steps`` controls boundary sampling only). If the
    fill looks blocky, increase ``pH_steps`` / ``Eh_steps``; the lines stay smooth
    but the colored regions resolve better.

    ``progress``
    ------------
    When ``True``, prints percent complete while scanning the pH–Eh grid.

    Junction labels
    ---------------
    Analytic junctions receive ``P1``, ``P2``, … in ``PourbaixResult.junction_points`` (dominant
    regions) and ``analytic_junction_points`` (all thermodynamic intersections).
    Use ``junction_table()`` or ``junction_coords("P4")`` for coordinates.
    ``plot(junction_label_style=...)`` controls annotations:

    - ``"numbered"`` — ``P1`` only on the figure; read coordinates from the result.
    - ``"numbered_coords"`` (default) — ``P1(7.00, 0.45)`` on the figure.
    - ``"species"`` — species-based label (``HSeO3- · SeO3-2 · …``).

    Boundary drawing
    ----------------
    ``plot(boundary_mode="dominant")`` (default) draws smooth analytic equilibrium curves
    only where two species share a border on the predominance diagram — not every place
    those species are equal in concentration. Junction markers use triple points of those
    regions only.

    ``plot(boundary_mode="all")`` draws every inferred analytic boundary and all line
    intersections (previous behaviour).

    Speciation methods
    ------------------
    ``model`` (default)
        Graph-based dominance using species-resolved Nernst boundaries, acid-base
        speciation, optional oligomer/hydration within levels, and Ksp post-checks.
        Fast and intended for exam-style diagrams at fixed total concentration.
        Saves analytic boundary coordinates and auto-labeled junction points.

        **Use when:** connected redox ladder per element, known E°/K/pKa, fixed
        ``element_totals``, and you need speed or exported line/point coordinates.

        **Limits:** ideal dilute Nernst + pKa; independent per-element chains; no
        cross-element redox; branched Latimer trees need manual wiring; oligomer
        regions depend on ``element_totals``; partial coupling for env16-style nets.

    ``equilibrium``
        Calls ``env.equilibrium()`` at every grid point with graph-based warm start.
        Primary saved geometry is the dominance grid matrix (``geometry_source='grid'``).

        **Use when:** graph model is incomplete or full reaction coupling is required.

        **Limits:** much slower; may fail on stiff networks.
    """

    environment: object
    track_species: Optional[Sequence[str]] = None
    element_totals: Optional[dict[str, float]] = None
    background_ions: Optional[dict[str, float]] = None
    junction_labels: Optional[dict[tuple[str, ...], str]] = None
    junction_label_style: JunctionLabelStyle = "numbered_coords"
    speciation_method: Literal["model", "equilibrium"] = "model"
    pH_min: float = 0.0
    pH_max: float = 14.0
    pH_steps: int = 40
    Eh_min: float = -1.0
    Eh_max: float = 1.5
    Eh_steps: int = 40
    geometry_pH_steps: int = 200
    equilibrium_method: str = "bgd"
    equilibrium_tol: float = 1e-8
    quotient_error_limit: float = 0.05
    min_concentration: float = 1e-30
    max_iter: int = 3000
    learning_rate: float = 0.3
    progress: bool = False

    def run(self) -> PourbaixResult:
        if _find_h_plus_index(self.environment) is None:
            raise ValueError("Pourbaix requires H+ in the environment.")

        graph = build_pourbaix_graph(self.environment)
        track_species = list(self.track_species or graph.all_track_species)
        if not track_species:
            raise ValueError("No trackable species found for Pourbaix.")

        totals = dict(self.element_totals or element_totals_from_env(self.environment, graph))

        pH_values = _pH_grid(self.pH_min, self.pH_max, self.pH_steps)
        eh_values = _eh_grid(self.Eh_min, self.Eh_max, self.Eh_steps)
        n_pH = len(pH_values)
        n_eh = len(eh_values)
        n_track = len(track_species)

        half_reactions = sorted(
            getattr(self.environment, "half_reactions", None) or [],
            key=lambda hr: hr.E0_SHE,
            reverse=True,
        )
        boundary_lines = [hr.boundary_line(pH_values) for hr in half_reactions]

        dominant = np.zeros((n_pH, n_eh), dtype=int)
        speciation = np.zeros((n_pH, n_eh, n_track), dtype=float)
        converged = np.ones((n_pH, n_eh), dtype=bool)

        label_to_index = {name: index for index, name in enumerate(track_species)}

        total_points = n_pH * n_eh
        done = 0
        for i, pH in enumerate(pH_values):
            for j, eh in enumerate(eh_values):
                if self.speciation_method == "equilibrium":
                    env_copy = self.environment.copy()
                    overrides: dict[str, float] = {}
                    neighbor = _neighbor_concentrations(speciation, track_species, i, j)
                    if neighbor is not None:
                        overrides.update(neighbor)
                    else:
                        guess_conc, _ = graph_speciation(
                            float(pH),
                            float(eh),
                            graph,
                            totals,
                            background_ions=self.background_ions,
                        )
                        overrides = dict(guess_conc)
                    apply_pourbaix_grid_state(env_copy, pH, eh, overrides)
                    try:
                        result = env_copy.equilibrium(
                            method=self.equilibrium_method,
                            tol=self.equilibrium_tol,
                            quotient_error_limit=self.quotient_error_limit,
                            min_concentration=self.min_concentration,
                            max_iter=self.max_iter,
                            learning_rate=self.learning_rate,
                            return_details=True,
                        )
                        concentrations = result.concentrations
                        converged[i, j] = bool(result.converged)
                        dominant[i, j] = _dominant_index(concentrations, env_copy, track_species)
                        label_map = {formula: idx for idx, formula in enumerate(env_copy.compound_labels)}
                        for k, label in enumerate(track_species):
                            idx = label_map.get(label)
                            if idx is not None:
                                speciation[i, j, k] = float(concentrations[idx])
                    except Exception:
                        converged[i, j] = False
                        dominant[i, j] = 0
                else:
                    concentrations, dom_index = graph_speciation(
                        float(pH),
                        float(eh),
                        graph,
                        totals,
                        background_ions=self.background_ions,
                    )
                    dom_name = graph.all_track_species[dom_index]
                    dominant[i, j] = label_to_index.get(dom_name, dom_index if dom_index < n_track else 0)
                    for k, label in enumerate(track_species):
                        speciation[i, j, k] = concentrations.get(label, 0.0)

                done += 1
                if self.progress and done % max(total_points // 10, 1) == 0:
                    print(f"Pourbaix grid: {100 * done / total_points:.0f}%")

        analytic_boundaries: list[PourbaixBoundary] = []
        analytic_junction_points: list[PourbaixJunction] = []
        geometry_source: Literal["analytic", "grid"] = "grid"

        if self.speciation_method == "model":
            analytic_boundaries, analytic_junction_points = compute_analytic_geometry(
                graph,
                pH_min=self.pH_min,
                pH_max=self.pH_max,
                pH_steps=self.geometry_pH_steps,
                eh_min=self.Eh_min,
                eh_max=self.Eh_max,
                junction_labels=self.junction_labels,
            )
            geometry_source = "analytic"

        equal_lines = _equal_concentration_boundaries(
            pH_values,
            eh_values,
            speciation,
            dominant,
            track_species,
            analytic_boundaries=analytic_boundaries or None,
            graph=graph,
            totals=totals,
            pH_min=self.pH_min,
            pH_max=self.pH_max,
            eh_min=self.Eh_min,
            eh_max=self.Eh_max,
            background_ions=self.background_ions,
        )
        dominant_junctions = _dominant_junction_points(
            pH_values,
            eh_values,
            dominant,
            track_species,
            junction_labels=self.junction_labels,
        )

        return PourbaixResult(
            grid_pH=pH_values,
            grid_Eh=eh_values,
            grid_dominant=dominant,
            track_species=track_species,
            speciation=speciation,
            converged=converged,
            geometry_source=geometry_source,
            analytic_boundaries=analytic_boundaries,
            junction_points=dominant_junctions,
            analytic_junction_points=analytic_junction_points,
            graph=graph,
            boundary_lines=boundary_lines,
            equal_boundary_lines=equal_lines,
        )
