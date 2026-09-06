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
PlotStyle = Literal["filled", "labeled"]
FrameEdge = Literal["pH_min", "pH_max", "Eh_min", "Eh_max"]


@dataclass(frozen=True)
class PourbaixFrameIntersection:
    pH: float
    Eh: float
    edge: FrameEdge


def _point_on_frame(
    pH: float,
    Eh: float,
    pH_min: float,
    pH_max: float,
    eh_min: float,
    eh_max: float,
    *,
    tol: float = 1e-9,
) -> Optional[FrameEdge]:
    on_ph_min = abs(pH - pH_min) <= tol and eh_min - tol <= Eh <= eh_max + tol
    on_ph_max = abs(pH - pH_max) <= tol and eh_min - tol <= Eh <= eh_max + tol
    on_eh_min = abs(Eh - eh_min) <= tol and pH_min - tol <= pH <= pH_max + tol
    on_eh_max = abs(Eh - eh_max) <= tol and pH_min - tol <= pH <= pH_max + tol
    hits = [name for flag, name in (
        (on_ph_min, "pH_min"),
        (on_ph_max, "pH_max"),
        (on_eh_min, "Eh_min"),
        (on_eh_max, "Eh_max"),
    ) if flag]
    if not hits:
        return None
    return hits[0]  # type: ignore[return-value]


def _segment_frame_intersections(
    pH0: float,
    eh0: float,
    pH1: float,
    eh1: float,
    pH_min: float,
    pH_max: float,
    eh_min: float,
    eh_max: float,
) -> list[PourbaixFrameIntersection]:
    points: list[PourbaixFrameIntersection] = []
    dpH = pH1 - pH0
    dEh = eh1 - eh0

    for edge, value, axis in (
        ("pH_min", pH_min, "pH"),
        ("pH_max", pH_max, "pH"),
        ("Eh_min", eh_min, "Eh"),
        ("Eh_max", eh_max, "Eh"),
    ):
        if axis == "pH":
            if abs(dpH) < 1e-15:
                if abs(pH0 - value) <= 1e-9:
                    for eh in (eh0, eh1):
                        if eh_min - 1e-9 <= eh <= eh_max + 1e-9:
                            points.append(PourbaixFrameIntersection(pH=value, Eh=float(eh), edge=edge))  # type: ignore[arg-type]
                continue
            t = (value - pH0) / dpH
            if -1e-9 <= t <= 1.0 + 1e-9:
                t = float(np.clip(t, 0.0, 1.0))
                eh = eh0 + t * dEh
                if eh_min - 1e-9 <= eh <= eh_max + 1e-9:
                    points.append(PourbaixFrameIntersection(pH=value, Eh=float(eh), edge=edge))  # type: ignore[arg-type]
        else:
            if abs(dEh) < 1e-15:
                if abs(eh0 - value) <= 1e-9:
                    for pH in (pH0, pH1):
                        if pH_min - 1e-9 <= pH <= pH_max + 1e-9:
                            points.append(PourbaixFrameIntersection(pH=float(pH), Eh=value, edge=edge))  # type: ignore[arg-type]
                continue
            t = (value - eh0) / dEh
            if -1e-9 <= t <= 1.0 + 1e-9:
                t = float(np.clip(t, 0.0, 1.0))
                pH = pH0 + t * dpH
                if pH_min - 1e-9 <= pH <= pH_max + 1e-9:
                    points.append(PourbaixFrameIntersection(pH=float(pH), Eh=value, edge=edge))  # type: ignore[arg-type]
    return points


def _polyline_frame_intersections(
    pH_line: np.ndarray,
    eh_line: np.ndarray,
    pH_min: float,
    pH_max: float,
    eh_min: float,
    eh_max: float,
) -> list[PourbaixFrameIntersection]:
    if len(pH_line) == 0:
        return []
    points: list[PourbaixFrameIntersection] = []
    for index in range(len(pH_line)):
        pH = float(pH_line[index])
        eh = float(eh_line[index])
        edge = _point_on_frame(pH, eh, pH_min, pH_max, eh_min, eh_max)
        if edge is not None:
            points.append(PourbaixFrameIntersection(pH=pH, Eh=eh, edge=edge))
    for index in range(len(pH_line) - 1):
        points.extend(
            _segment_frame_intersections(
                float(pH_line[index]),
                float(eh_line[index]),
                float(pH_line[index + 1]),
                float(eh_line[index + 1]),
                pH_min,
                pH_max,
                eh_min,
                eh_max,
            )
        )
    return _dedupe_frame_intersections(points)


def _dedupe_frame_intersections(
    points: Sequence[PourbaixFrameIntersection],
    *,
    pH_decimals: int = 3,
    eh_decimals: int = 3,
) -> list[PourbaixFrameIntersection]:
    seen: set[tuple[float, float]] = set()
    unique: list[PourbaixFrameIntersection] = []
    for point in points:
        key = (round(point.pH, pH_decimals), round(point.Eh, eh_decimals))
        if key in seen:
            continue
        seen.add(key)
        unique.append(point)
    return sorted(unique, key=lambda item: (item.pH, item.Eh))


def _plotted_boundary_polylines(
    result: "PourbaixResult",
    *,
    boundary_mode: BoundaryMode,
    show_equal_boundaries: bool,
    show_analytic_boundaries: bool,
    include_water_lines: bool,
) -> list[tuple[np.ndarray, np.ndarray]]:
    polylines: list[tuple[np.ndarray, np.ndarray]] = []
    if show_equal_boundaries and result.equal_boundary_lines:
        for _, _, pH_line, eh_line in result.equal_boundary_lines:
            polylines.append((pH_line, eh_line))
    if show_analytic_boundaries and result.analytic_boundaries:
        for boundary in result.analytic_boundaries:
            if boundary.kind == "water" and not include_water_lines:
                continue
            polylines.append((boundary.pH, boundary.Eh))
    if include_water_lines and not (
        show_analytic_boundaries and any(b.kind == "water" for b in result.analytic_boundaries)
    ):
        polylines.append((result.grid_pH, 1.229 - 2 * NERNST_K * result.grid_pH))
        polylines.append((result.grid_pH, -2 * NERNST_K * result.grid_pH))
    return polylines


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
    grid_pH: Optional[np.ndarray] = None,
) -> bool:
    """True when (pH, Eh) lies on the border between two dominant region species."""
    label_to_index = {name: index for index, name in enumerate(track_species)}

    if grid_pH is not None and len(grid_pH) > 1:
        pH_delta = max(0.05, 0.45 * (float(grid_pH[-1]) - float(grid_pH[0])) / (len(grid_pH) - 1))
    else:
        pH_delta = 0.05
    eh_delta = 1e-3

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
        dominant_index(pH, Eh + eh_delta),
        dominant_index(pH, Eh - eh_delta),
        dominant_index(pH + pH_delta, Eh),
        dominant_index(pH - pH_delta, Eh),
    }
    seen.discard(-1)
    return idx_a in seen and idx_b in seen


def _merge_eh_intervals(ehs: Sequence[float], eh_gap: float) -> list[tuple[float, float]]:
    if not ehs:
        return []
    ordered = sorted(set(float(eh) for eh in ehs))
    intervals: list[tuple[float, float]] = []
    start = ordered[0]
    previous = ordered[0]
    for eh in ordered[1:]:
        if eh - previous > eh_gap:
            intervals.append((start, previous))
            start = eh
        previous = eh
    intervals.append((start, previous))
    return intervals


def _vertical_eh_spans_for_pair(
    grid_pH: np.ndarray,
    grid_Eh: np.ndarray,
    grid_dominant: np.ndarray,
    idx_a: int,
    idx_b: int,
    *,
    pH_target: Optional[float] = None,
    pH_tol: Optional[float] = None,
) -> list[tuple[float, float, float]]:
    """Return ``(pH_line, eh_lo, eh_hi)`` segments for a vertical predominance border."""
    if pH_tol is None:
        pH_tol = (float(grid_pH[-1]) - float(grid_pH[0])) / max(len(grid_pH) - 1, 1) * 0.75
    eh_gap = (float(grid_Eh[-1]) - float(grid_Eh[0])) / max(len(grid_Eh) - 1, 1) * 1.5

    edges: dict[float, list[float]] = {}
    for i in range(len(grid_pH) - 1):
        pH_edge = 0.5 * (float(grid_pH[i]) + float(grid_pH[i + 1]))
        if pH_target is not None and abs(pH_edge - pH_target) > pH_tol:
            continue
        for j in range(len(grid_Eh)):
            left = int(grid_dominant[i, j])
            right = int(grid_dominant[i + 1, j])
            if {left, right} != {idx_a, idx_b}:
                continue
            key = round(pH_edge, 4)
            edges.setdefault(key, []).append(float(grid_Eh[j]))

    segments: list[tuple[float, float, float]] = []
    for pH_edge, ehs in edges.items():
        pH_line = float(pH_target) if pH_target is not None else pH_edge
        for eh_lo, eh_hi in _merge_eh_intervals(ehs, eh_gap):
            if eh_hi >= eh_lo:
                segments.append((pH_line, eh_lo, eh_hi))
    return segments


def _dominant_acid_base_lines(
    grid_pH: np.ndarray,
    grid_Eh: np.ndarray,
    grid_dominant: np.ndarray,
    track_species: Sequence[str],
    analytic_boundaries: Sequence[PourbaixBoundary],
    *,
    pH_min: float,
    pH_max: float,
    eh_min: float,
    eh_max: float,
) -> list[tuple[int, int, np.ndarray, np.ndarray]]:
    """Vertical pH boundaries between dominant acid/base neighbors."""
    lines: list[tuple[int, int, np.ndarray, np.ndarray]] = []
    acid_base_targets = {
        frozenset({boundary.left, boundary.right})
        for boundary in analytic_boundaries
        if boundary.kind == "acid_base"
    }
    use_grid_only = not acid_base_targets

    for idx_a, idx_b in sorted(_adjacent_species_pairs(grid_dominant)):
        name_a = track_species[idx_a]
        name_b = track_species[idx_b]
        target = frozenset({name_a, name_b})
        if not use_grid_only and target not in acid_base_targets:
            continue

        pH_target: Optional[float] = None
        if not use_grid_only:
            for boundary in analytic_boundaries:
                if boundary.kind == "acid_base" and {boundary.left, boundary.right} == {name_a, name_b}:
                    pH_target = float(boundary.pH[0])
                    break

        spans = _vertical_eh_spans_for_pair(
            grid_pH,
            grid_Eh,
            grid_dominant,
            idx_a,
            idx_b,
            pH_target=pH_target,
        )
        if use_grid_only and not spans:
            continue

        for pH_line, eh_lo, eh_hi in spans:
            eh_lo = max(eh_lo, eh_min)
            eh_hi = min(eh_hi, eh_max)
            if eh_hi <= eh_lo:
                continue
            lines.append(
                (
                    idx_a,
                    idx_b,
                    np.asarray([pH_line, pH_line], dtype=float),
                    np.asarray([eh_lo, eh_hi], dtype=float),
                )
            )
    return lines


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
    grid_pH: Optional[np.ndarray] = None,
) -> list[tuple[np.ndarray, np.ndarray]]:
    if boundary.kind == "acid_base":
        return []
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
            grid_pH=grid_pH,
        ):
            continue
        active[index] = True
    return _split_active_segments(pH_curve, eh_curve, active)


def _dominant_analytic_boundary_lines(
    grid_pH: np.ndarray,
    grid_Eh: np.ndarray,
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
            if boundary.kind in {"water", "acid_base"}:
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
                grid_pH=grid_pH,
            ):
                if len(pH_seg) >= 2:
                    lines.append((idx_a, idx_b, pH_seg, eh_seg))
    lines.extend(
        _dominant_acid_base_lines(
            grid_pH,
            grid_Eh,
            grid_dominant,
            track_species,
            analytic_boundaries,
            pH_min=pH_min,
            pH_max=pH_max,
            eh_min=eh_min,
            eh_max=eh_max,
        )
    )
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
                                grid_pH=grid_pH,
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
            grid_pH,
            grid_Eh,
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


def _connected_dominant_regions(
    grid_pH: np.ndarray,
    grid_Eh: np.ndarray,
    grid_dominant: np.ndarray,
) -> list[tuple[int, list[tuple[int, int]]]]:
    """Four-connected components of equal ``grid_dominant`` index."""
    n_pH, n_eh = grid_dominant.shape
    visited = np.zeros((n_pH, n_eh), dtype=bool)
    components: list[tuple[int, list[tuple[int, int]]]] = []
    for i in range(n_pH):
        for j in range(n_eh):
            if visited[i, j]:
                continue
            species = int(grid_dominant[i, j])
            cells: list[tuple[int, int]] = []
            stack = [(i, j)]
            visited[i, j] = True
            while stack:
                ci, cj = stack.pop()
                cells.append((ci, cj))
                for ni, nj in ((ci - 1, cj), (ci + 1, cj), (ci, cj - 1), (ci, cj + 1)):
                    if 0 <= ni < n_pH and 0 <= nj < n_eh and not visited[ni, nj]:
                        if int(grid_dominant[ni, nj]) != species:
                            continue
                        visited[ni, nj] = True
                        stack.append((ni, nj))
            components.append((species, cells))
    return components


def _region_label_fontsize(
    ax,
    text: str,
    pH_lo: float,
    pH_hi: float,
    eh_lo: float,
    eh_hi: float,
    *,
    min_font: float = 6.0,
    max_font: float = 14.0,
) -> float:
    """Pick a font size that fits inside a predominance region bbox."""
    fig = ax.figure
    dpi = fig.dpi
    transform = ax.transData.transform
    (x0, y0) = transform((pH_lo, eh_lo))
    (x1, y1) = transform((pH_hi, eh_hi))
    width_pt = abs(x1 - x0) * 72.0 / dpi
    height_pt = abs(y1 - y0) * 72.0 / dpi
    if width_pt <= 0 or height_pt <= 0:
        return min_font
    char_count = max(len(text), 1)
    by_width = width_pt / (0.55 * char_count)
    by_height = height_pt * 0.45
    return float(np.clip(min(by_width, by_height), min_font, max_font))


def _dominant_region_annotations(
    grid_pH: np.ndarray,
    grid_Eh: np.ndarray,
    grid_dominant: np.ndarray,
    label_map: dict[int, str],
    *,
    ax,
    min_cells: int = 2,
) -> None:
    """Place species labels at connected-region centroids with size-aware font."""
    if len(grid_pH) < 2 or len(grid_Eh) < 2:
        return
    dpH = float(grid_pH[1] - grid_pH[0])
    dEh = float(grid_Eh[1] - grid_Eh[0])
    for species_index, cells in _connected_dominant_regions(grid_pH, grid_Eh, grid_dominant):
        if len(cells) < min_cells:
            continue
        text = label_map.get(species_index)
        if not text:
            continue
        pH_values = [float(grid_pH[i]) for i, _ in cells]
        eh_values = [float(grid_Eh[j]) for _, j in cells]
        pH_center = float(np.mean(pH_values))
        eh_center = float(np.mean(eh_values))
        pH_lo = min(pH_values) - 0.5 * dpH
        pH_hi = max(pH_values) + 0.5 * dpH
        eh_lo = min(eh_values) - 0.5 * dEh
        eh_hi = max(eh_values) + 0.5 * dEh
        fontsize = _region_label_fontsize(ax, text, pH_lo, pH_hi, eh_lo, eh_hi)
        ax.text(
            pH_center,
            eh_center,
            text,
            ha="center",
            va="center",
            fontsize=fontsize,
            clip_on=True,
            wrap=True,
        )


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

    def frame_intersections(
        self,
        *,
        boundary_mode: BoundaryMode = "dominant",
        show_equal_boundaries: Optional[bool] = None,
        show_analytic_boundaries: Optional[bool] = None,
        include_water_lines: bool = True,
    ) -> list[PourbaixFrameIntersection]:
        """Points where plotted boundaries meet the diagram frame (min/max pH and Eh)."""
        if show_analytic_boundaries is None:
            show_analytic_boundaries = boundary_mode == "all"
        if show_equal_boundaries is None:
            show_equal_boundaries = boundary_mode == "dominant"

        pH_min = float(self.grid_pH[0])
        pH_max = float(self.grid_pH[-1])
        eh_min = float(self.grid_Eh[0])
        eh_max = float(self.grid_Eh[-1])
        points: list[PourbaixFrameIntersection] = []
        for pH_line, eh_line in _plotted_boundary_polylines(
            self,
            boundary_mode=boundary_mode,
            show_equal_boundaries=show_equal_boundaries,
            show_analytic_boundaries=show_analytic_boundaries,
            include_water_lines=include_water_lines,
        ):
            points.extend(
                _polyline_frame_intersections(
                    np.asarray(pH_line, dtype=float),
                    np.asarray(eh_line, dtype=float),
                    pH_min,
                    pH_max,
                    eh_min,
                    eh_max,
                )
            )
        return _dedupe_frame_intersections(points)

    def frame_intersection_table(
        self,
        *,
        boundary_mode: BoundaryMode = "dominant",
        show_equal_boundaries: Optional[bool] = None,
        show_analytic_boundaries: Optional[bool] = None,
        include_water_lines: bool = True,
    ) -> list[dict]:
        return [
            {"pH": point.pH, "Eh": point.Eh, "edge": point.edge}
            for point in self.frame_intersections(
                boundary_mode=boundary_mode,
                show_equal_boundaries=show_equal_boundaries,
                show_analytic_boundaries=show_analytic_boundaries,
                include_water_lines=include_water_lines,
            )
        ]

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
        plot_style: PlotStyle = "filled",
        boundary_mode: BoundaryMode = "dominant",
        show_analytic_boundaries: Optional[bool] = None,
        show_equal_boundaries: Optional[bool] = None,
        show_junction_labels: bool = True,
        junction_label_style: JunctionLabelStyle = "numbered_coords",
        junction_pH_decimals: int = 2,
        junction_eh_decimals: int = 2,
        show_frame_intersections: bool = False,
        frame_intersection_labels: bool = True,
        frame_intersection_pH_decimals: int = 2,
        frame_intersection_eh_decimals: int = 2,
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

        label_map = labels or {i: name for i, name in enumerate(self.track_species)}

        if plot_style == "filled":
            pH_mesh, eh_mesh = np.meshgrid(self.grid_pH, self.grid_Eh, indexing="ij")
            n_species = max(len(self.track_species), int(self.grid_dominant.max()) + 1)
            cmap = ListedColormap(plt.cm.tab10.colors[: max(n_species, 1)])
            norm = BoundaryNorm(np.arange(-0.5, n_species + 0.5, 1), n_species)
            ax.pcolormesh(pH_mesh, eh_mesh, self.grid_dominant, cmap=cmap, norm=norm, shading="auto")
        else:
            ax.set_facecolor("white")

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

        if show_frame_intersections:
            for point in self.frame_intersections(
                boundary_mode=boundary_mode,
                show_equal_boundaries=show_equal_boundaries,
                show_analytic_boundaries=show_analytic_boundaries,
                include_water_lines=include_water_lines,
            ):
                ax.plot(point.pH, point.Eh, "ks", ms=5, mfc="none")
                if frame_intersection_labels:
                    ax.annotate(
                        f"({point.pH:.{frame_intersection_pH_decimals}f}, "
                        f"{point.Eh:.{frame_intersection_eh_decimals}f})",
                        (point.pH, point.Eh),
                        fontsize=5,
                        xytext=(3, -8),
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

        if plot_style == "labeled":
            _dominant_region_annotations(
                self.grid_pH,
                self.grid_Eh,
                self.grid_dominant,
                label_map,
                ax=ax,
            )
        else:
            n_species = max(len(self.track_species), int(self.grid_dominant.max()) + 1)
            cmap = ListedColormap(plt.cm.tab10.colors[: max(n_species, 1)])
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
        show_frame_intersections: bool = False,
        frame_intersection_labels: bool = True,
        frame_intersection_pH_decimals: int = 2,
        frame_intersection_eh_decimals: int = 2,
    ):
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(8, 6))
        show_equal = boundary_mode == "dominant"
        show_analytic = boundary_mode == "all"
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
        if show_frame_intersections:
            for point in self.frame_intersections(
                boundary_mode=boundary_mode,
                show_equal_boundaries=show_equal,
                show_analytic_boundaries=show_analytic,
                include_water_lines=include_water_lines,
            ):
                ax.plot(point.pH, point.Eh, "ks", ms=5, mfc="none")
                if frame_intersection_labels:
                    ax.annotate(
                        f"({point.pH:.{frame_intersection_pH_decimals}f}, "
                        f"{point.Eh:.{frame_intersection_eh_decimals}f})",
                        (point.pH, point.Eh),
                        fontsize=6,
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
    ``pH_steps`` and ``Eh_steps`` set the **predominance-region grid** used to decide
    which species border each other and (for ``plot_style="filled"``) the colored
    background. Each cell stores one dominant species index.

    +------------------+---------------------------+---------------------------+
    | Parameter        | ``speciation_method``     | Primary effect            |
    +==================+===========================+===========================+
    | ``pH_steps``     | both                      | Region layout + labels    |
    | ``Eh_steps``     | both                      | Region layout + labels    |
    | ``geometry_pH_steps`` | ``model`` only       | Analytic boundary sampling|
    +------------------+---------------------------+---------------------------+

    **``model``:** boundaries are analytic Nernst/pKa curves; ``geometry_pH_steps``
    controls line smoothness. ``pH_steps`` / ``Eh_steps`` mainly affect which
    neighbors are detected and how blocky the fill is — use ~25–40 for publication
    fill, or coarser for ``plot_style="labeled"``.

    **``equilibrium``:** no analytic geometry; ``pH_steps`` / ``Eh_steps`` control
    both region assignment and boundary placement (finer = more accurate, slower).

    Plot styles
    -----------
    ``plot_style="filled"`` (default) — colored regions plus boundary lines.

    ``plot_style="labeled"`` — white background, boundary lines, and each connected
    predominant region annotated with its species at the region centroid (font size
    scales with region area). No color legend.

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

    Frame intersections
    -------------------
    ``plot(show_frame_intersections=True)`` marks open squares where drawn boundaries
    meet the diagram frame (minimum/maximum pH and Eh). Query coordinates with
    ``frame_intersection_table()`` on the result object.

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
