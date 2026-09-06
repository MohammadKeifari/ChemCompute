"""Graph-based Pourbaix speciation, boundaries, and analytic geometry."""

from __future__ import annotations

import math
import re
import warnings
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Literal, Optional, Sequence

import numpy as np

from ._general import POURBAIX_RESERVED_SPECIES
from ._half_reaction import redox_reaction_indices

NERNST_K = 0.05916
POURBAIX_SKIP_SPECIES = POURBAIX_RESERVED_SPECIES | frozenset({"H2O"})


def _formula(entry) -> str:
    compound = entry["compound"]
    return compound.formula if hasattr(compound, "formula") else str(compound)


def _primary_side_species(entries, *, skip: frozenset[str]) -> str:
    for entry in entries:
        formula = _formula(entry)
        if formula not in skip:
            return formula
    raise ValueError("Could not identify a primary species for Pourbaix inference.")


def _build_species_chain(start: str, pka_pairs: Sequence[tuple[str, str, float]]) -> list[str]:
    """Return the acid→base path that contains ``start`` (most acidic first)."""
    acid_to_base = {acid: base for acid, base, _ in pka_pairs}
    base_to_acid = {base: acid for acid, base, _ in pka_pairs}
    current = start
    seen = {current}
    while current in base_to_acid:
        acid = base_to_acid[current]
        if acid in seen:
            break
        current = acid
        seen.add(current)
    chain = [current]
    while current in acid_to_base:
        nxt = acid_to_base[current]
        if nxt in chain:
            break
        chain.append(nxt)
        current = nxt
    return chain


def _extract_pka_pairs(env) -> list[tuple[str, str, float]]:
    redox = set(redox_reaction_indices(env))
    pairs: list[tuple[str, str, float]] = []
    for index, rxn in enumerate(env.reactions):
        if index in redox:
            continue
        reactants = rxn.reactants
        products = rxn.products
        if len(reactants) != 1 or len(products) != 2:
            continue
        acid = _formula(reactants[0])
        if acid in POURBAIX_SKIP_SPECIES:
            continue
        product_formulas = {_formula(entry): entry for entry in products}
        if "H+" not in product_formulas:
            continue
        base_candidates = [name for name in product_formulas if name not in POURBAIX_SKIP_SPECIES]
        if len(base_candidates) != 1:
            continue
        base = base_candidates[0]
        pka = -math.log10(max(float(rxn.K), 1e-300))
        pairs.append((acid, base, pka))
    return pairs


def _acid_base_form(pH: float, candidates: Sequence[str], pka_pairs: Sequence[tuple[str, str, float]]) -> str:
    active = candidates[0]
    candidate_set = set(candidates)
    for acid, base, pka in sorted(pka_pairs, key=lambda item: item[2]):
        if acid in candidate_set and base in candidate_set and pH >= pka:
            active = base
    return active


def _chain_concentrations(
    pH: float,
    species: Sequence[str],
    pka_pairs: Sequence[tuple[str, str, float]],
    total: float,
) -> dict[str, float]:
    if total <= 0:
        return {name: 0.0 for name in species}
    if len(species) == 1:
        return {species[0]: total}
    h = 10.0 ** (-float(pH))
    pair_map = {(acid, base): pka for acid, base, pka in pka_pairs}
    weights = [1.0]
    for index in range(len(species) - 1):
        key = (species[index], species[index + 1])
        pka = pair_map.get(key)
        if pka is None:
            weights.append(weights[-1])
            continue
        weights.append(weights[-1] * (10 ** (-pka) / max(h, 1e-300)))
    norm = sum(weights)
    return {species[index]: total * weights[index] / norm for index in range(len(species))}


def _boundary_e(segments: Sequence[tuple[float, float, float, float]], pH: float) -> float:
    pH = float(pH)
    for ph_hi, e0, m, n in segments:
        if pH <= ph_hi:
            return e0 - NERNST_K * m / n * pH
    ph_hi, e0, m, n = segments[-1]
    return e0 - NERNST_K * m / n * pH


def _build_segments_reduced(hr, species_chain, pka_pairs):
    pair_map = {(acid, base): pka for acid, base, pka in pka_pairs}
    n = hr.n_electrons
    m = hr.net_h_plus_stoichiometry()
    e0 = hr.E0_SHE
    segments: list[tuple[float, float, float, float]] = []
    for index in range(len(species_chain) - 1):
        acid = species_chain[index]
        base = species_chain[index + 1]
        pka = pair_map.get((acid, base))
        if pka is None:
            continue
        e_at = e0 - NERNST_K * m / n * pka
        segments.append((pka, e0, m, n))
        m = max(m - 1.0, 0.0)
        e0 = e_at + NERNST_K * m / n * pka
    segments.append((14.0, e0, m, n))
    return segments


def _build_segments(hr, species_chain, pka_pairs):
    pair_map = {(acid, base): pka for acid, base, pka in pka_pairs}
    n = hr.n_electrons
    m = hr.net_h_plus_stoichiometry()
    e0 = hr.E0_SHE
    segments: list[tuple[float, float, float, float]] = []
    for index in range(len(species_chain) - 1):
        acid = species_chain[index]
        base = species_chain[index + 1]
        pka = pair_map.get((acid, base))
        if pka is None:
            continue
        e_at = e0 - NERNST_K * m / n * pka
        segments.append((pka, e0, m, n))
        m += 1.0
        e0 = e_at + NERNST_K * m / n * pka
    segments.append((14.0, e0, m, n))
    return segments

BoundaryKind = Literal["redox", "acid_base", "ksp", "association", "water"]

_SKIP_ELEMENTS = frozenset({"H", "O"})
_ELEMENT_RE = re.compile(r"[A-Z][a-z]?")


def _redox_element(formula: str) -> str:
    """Return the redox-active element symbol for a species formula."""
    base = formula.split(".")[0]
    for match in _ELEMENT_RE.finditer(base):
        symbol = match.group(0)
        if symbol not in _SKIP_ELEMENTS:
            return symbol
    raise ValueError(f"Could not identify redox element in formula {formula!r}.")


def _count_element(formula: str, element: str) -> int:
    """Count atoms of ``element`` in a species formula (charge/phase stripped)."""
    core = formula.split(".")[0]
    core = re.sub(r"[+-]\d*$", "", core)
    total = 0
    for match in _ELEMENT_RE.finditer(core):
        if match.group(0) != element:
            continue
        digits = re.match(r"\d+", core[match.end() :])
        total += int(digits.group(0)) if digits else 1
    return total


def _primary_with_coeff(entries) -> tuple[str, float]:
    for entry in entries:
        formula = _formula(entry)
        if formula not in POURBAIX_SKIP_SPECIES:
            return formula, float(entry.get("stoichiometric_coefficient", 1.0))
    raise ValueError("Could not identify a primary species for Pourbaix inference.")


def _edge_atom_balanced(hr, element: str) -> bool:
    ox, c_ox = _primary_with_coeff(hr.oxidized)
    red, c_red = _primary_with_coeff(hr.reduced)
    n_ox = _count_element(ox, element) * c_ox
    n_red = _count_element(red, element) * c_red
    return n_ox > 0 and abs(n_ox - n_red) < 1e-9


def format_junction_label(species: Sequence[str]) -> str:
    return " · ".join(sorted(species))


@dataclass
class AcidBaseEdge:
    acid: str
    base: str
    pka: float


@dataclass
class PrecipitationEdge:
    solid: str
    ions: dict[str, float]
    ksp: float
    reaction_index: int = -1


@dataclass
class AssociationEdge:
    reactants: dict[str, float]
    product: str
    k: float
    reaction_index: int = -1


@dataclass
class HydrationEdge:
    dehydrated: str
    hydrated: str
    kh: float
    reaction_index: int = -1


@dataclass
class ElementChain:
    element: str
    oxidation_levels: list[list[str]]
    half_reactions: list
    track_species: list[str]
    pka_pairs: list[tuple[str, str, float]]
    acid_base_edges: list[AcidBaseEdge] = field(default_factory=list)
    precipitation_edges: list[PrecipitationEdge] = field(default_factory=list)
    association_edges: list[AssociationEdge] = field(default_factory=list)
    hydration_edges: list[HydrationEdge] = field(default_factory=list)
    redox_edges: list[tuple[int, int]] = field(default_factory=list)


@dataclass
class PourbaixGraph:
    chains: list[ElementChain]
    all_track_species: list[str]
    precipitation_edges: list[PrecipitationEdge] = field(default_factory=list)
    acid_base_edges: list[AcidBaseEdge] = field(default_factory=list)
    association_edges: list[AssociationEdge] = field(default_factory=list)
    hydration_edges: list[HydrationEdge] = field(default_factory=list)


@dataclass
class PourbaixBoundary:
    left: str
    right: str
    kind: BoundaryKind
    pH: np.ndarray
    Eh: np.ndarray
    level_pair: tuple[int, int] = (-1, -1)
    chain_element: str = ""


@dataclass
class PourbaixJunction:
    pH: float
    Eh: float
    species: tuple[str, ...]
    label: str
    point_id: str = ""
    boundary_kinds: tuple[str, ...] = ()


JunctionLabelStyle = Literal["species", "numbered", "numbered_coords"]


def assign_junction_point_ids(junctions: list[PourbaixJunction]) -> list[PourbaixJunction]:
    """Sort junctions by pH then Eh and assign ``P1``, ``P2``, … to ``point_id``."""
    ordered = sorted(junctions, key=lambda item: (item.pH, item.Eh))
    for index, junction in enumerate(ordered, start=1):
        junction.point_id = f"P{index}"
    return ordered


def format_junction_plot_label(
    junction: PourbaixJunction,
    *,
    style: JunctionLabelStyle = "numbered_coords",
    pH_decimals: int = 2,
    eh_decimals: int = 2,
) -> str:
    """Build the annotation text for a junction marker on a plot."""
    point_id = junction.point_id or "P?"
    if style == "species":
        return junction.label
    if style == "numbered":
        return point_id
    return f"{point_id}({junction.pH:.{pH_decimals}f}, {junction.Eh:.{eh_decimals}f})"


def _is_solid(entry) -> bool:
    compound = entry["compound"]
    phase_list = getattr(compound, "phase_point_list", None) or []
    for point in phase_list:
        if point.get("phase") == "s":
            return True
    name = _formula(entry)
    return name.endswith(".s")


def _is_liquid_water(entry) -> bool:
    return _formula(entry) == "H2O"


def classify_reaction(rxn, index: int, *, redox_indices: set[int]) -> Optional[object]:
    """Classify a non-redox environment reaction into a graph edge."""
    if index in redox_indices:
        return None

    reactants = rxn.reactants
    products = rxn.products
    k = float(rxn.K)

    # Acid-base: HA ⇌ H+ + A-
    if len(reactants) == 1 and len(products) == 2:
        acid = _formula(reactants[0])
        if acid not in POURBAIX_SKIP_SPECIES:
            product_formulas = {_formula(entry): entry for entry in products}
            if "H+" in product_formulas:
                bases = [name for name in product_formulas if name not in POURBAIX_SKIP_SPECIES]
                if len(bases) == 1:
                    return AcidBaseEdge(acid, bases[0], -math.log10(max(k, 1e-300)))

    # Hydration: X + H2O ⇌ Y  or  Y ⇌ X + H2O
    for side_a, side_b, forward in ((reactants, products, True), (products, reactants, False)):
        if len(side_a) == 2 and len(side_b) == 1:
            formulas_a = [_formula(e) for e in side_a]
            if "H2O" in formulas_a or any(_is_liquid_water(e) for e in side_a):
                other = [name for name in formulas_a if name not in POURBAIX_SKIP_SPECIES and name != "H2O"]
                if len(other) == 1 and not _is_solid(side_b[0]):
                    dehydrated, hydrated = (other[0], _formula(side_b[0])) if forward else (_formula(side_b[0]), other[0])
                    return HydrationEdge(dehydrated, hydrated, k if forward else 1.0 / max(k, 1e-300), index)

    solids_r = [e for e in reactants if _is_solid(e)]
    solids_p = [e for e in products if _is_solid(e)]
    skip_ksp = frozenset({"H2O"})
    aq_r = [e for e in reactants if not _is_solid(e) and _formula(e) not in skip_ksp]
    aq_p = [e for e in products if not _is_solid(e) and _formula(e) not in skip_ksp]

    # Ksp: Solid ⇌ ions
    if len(solids_r) == 1 and not solids_p and aq_p:
        ions = {_formula(e): float(e["stoichiometric_coefficient"]) for e in aq_p}
        return PrecipitationEdge(_formula(solids_r[0]), ions, k, index)
    if len(solids_p) == 1 and not solids_r and aq_r:
        ions = {_formula(e): float(e["stoichiometric_coefficient"]) for e in aq_r}
        return PrecipitationEdge(_formula(solids_p[0]), ions, 1.0 / max(k, 1e-300), index)

    # Association: n A ⇌ A_n  or  A + B ⇌ AB
    if len(reactants) >= 2 and len(products) == 1 and not _is_solid(products[0]):
        rdict = {_formula(e): float(e["stoichiometric_coefficient"]) for e in reactants}
        if all(name not in POURBAIX_SKIP_SPECIES for name in rdict):
            return AssociationEdge(rdict, _formula(products[0]), k, index)
    if len(products) >= 2 and len(reactants) == 1 and not _is_solid(reactants[0]):
        pdict = {_formula(e): float(e["stoichiometric_coefficient"]) for e in products}
        if all(name not in POURBAIX_SKIP_SPECIES for name in pdict):
            return AssociationEdge(pdict, _formula(reactants[0]), 1.0 / max(k, 1e-300), index)

    return None


def _classify_all_reactions(env) -> tuple[list[AcidBaseEdge], list[PrecipitationEdge], list[AssociationEdge], list[HydrationEdge]]:
    redox = set(redox_reaction_indices(env))
    acid_base: list[AcidBaseEdge] = []
    precipitation: list[PrecipitationEdge] = []
    association: list[AssociationEdge] = []
    hydration: list[HydrationEdge] = []
    warned: set[int] = set()

    for index, rxn in enumerate(env.reactions):
        if index in redox:
            continue
        edge = classify_reaction(rxn, index, redox_indices=redox)
        if edge is None:
            continue
        if isinstance(edge, AcidBaseEdge):
            acid_base.append(edge)
        elif isinstance(edge, PrecipitationEdge):
            precipitation.append(edge)
        elif isinstance(edge, AssociationEdge):
            association.append(edge)
        elif isinstance(edge, HydrationEdge):
            hydration.append(edge)
        else:
            if index not in warned:
                warnings.warn(f"Reaction {index} not used in Pourbaix graph model.", stacklevel=2)
                warned.add(index)

    if not acid_base:
        acid_base = [AcidBaseEdge(a, b, p) for a, b, p in _extract_pka_pairs(env)]
    return acid_base, precipitation, association, hydration


def _merge_species_groups(groups: list[list[str]]) -> list[list[str]]:
    """Union groups that share any species, preserving first-seen order inside each group."""
    merged: list[list[str]] = []
    for group in groups:
        group_set = set(group)
        hits = [index for index, existing in enumerate(merged) if group_set & set(existing)]
        if not hits:
            merged.append(list(group))
            continue
        anchor = hits[0]
        for index in reversed(hits[1:]):
            for name in merged[index]:
                if name not in merged[anchor]:
                    merged[anchor].append(name)
            del merged[index]
        for name in group:
            if name not in merged[anchor]:
                merged[anchor].append(name)
    return merged


def _build_oxidation_graph(
    element: str,
    hrs: Sequence,
    pka_pairs: Sequence[tuple[str, str, float]],
    hydration: Sequence[HydrationEdge],
) -> tuple[list[list[str]], list, list[tuple[int, int]]]:
    """Group HR species into oxidation levels and return directed ox→red edges."""
    groups: list[list[str]] = []
    for hr in hrs:
        ox = _primary_side_species(hr.oxidized, skip=POURBAIX_SKIP_SPECIES)
        red = _primary_side_species(hr.reduced, skip=POURBAIX_SKIP_SPECIES)
        for name in (ox, red):
            groups.append(
                _extend_level_with_hydration(
                    _build_species_chain(name, pka_pairs),
                    pka_pairs,
                    hydration,
                )
            )
    levels = _merge_species_groups(groups)

    def level_of(name: str) -> int:
        for index, level in enumerate(levels):
            if name in level:
                return index
        raise ValueError(f"Species {name!r} is not in an oxidation level for {element}.")

    raw_edges: list[tuple[int, int, object]] = []
    seen_edge: set[tuple[int, int, int]] = set()
    for hr in hrs:
        ox = _primary_side_species(hr.oxidized, skip=POURBAIX_SKIP_SPECIES)
        red = _primary_side_species(hr.reduced, skip=POURBAIX_SKIP_SPECIES)
        ox_i = level_of(ox)
        red_i = level_of(red)
        if ox_i == red_i:
            continue
        key = (ox_i, red_i, id(hr))
        if key in seen_edge:
            continue
        seen_edge.add(key)
        raw_edges.append((ox_i, red_i, hr))

    reduced_names = {
        _primary_side_species(hr.reduced, skip=POURBAIX_SKIP_SPECIES) for hr in hrs
    }
    oxidized_names = {
        _primary_side_species(hr.oxidized, skip=POURBAIX_SKIP_SPECIES) for hr in hrs
    }
    roots = [
        index
        for index, level in enumerate(levels)
        if any(name in oxidized_names for name in level)
        and not any(name in reduced_names for name in level)
    ]
    if not roots:
        roots = [0]

    children: dict[int, list[tuple[float, int]]] = defaultdict(list)
    for ox_i, red_i, hr in raw_edges:
        red_name = _primary_side_species(hr.reduced, skip=POURBAIX_SKIP_SPECIES)
        n_atom = max(_count_element(red_name, element), 1)
        children[ox_i].append((hr.n_electrons / n_atom, red_i))

    ordered: list[int] = []
    seen_levels: set[int] = set()
    queue: deque[int] = deque(sorted(roots))
    while queue:
        current = queue.popleft()
        if current in seen_levels:
            continue
        seen_levels.add(current)
        ordered.append(current)
        for _, child in sorted(children.get(current, [])):
            if child not in seen_levels:
                queue.append(child)
    for index in range(len(levels)):
        if index not in seen_levels:
            ordered.append(index)

    old_to_new = {old: new for new, old in enumerate(ordered)}
    oxidation_levels = [levels[index] for index in ordered]

    tree_pairs: set[tuple[int, int]] = set()
    for parent in ordered:
        for _, child in sorted(children.get(parent, [])):
            tree_pairs.add((old_to_new[parent], old_to_new[child]))

    bfs_hrs: list = []
    bfs_edges: list[tuple[int, int]] = []
    leftover_hrs: list = []
    leftover_edges: list[tuple[int, int]] = []
    used_hr: set[int] = set()
    for ox_i, red_i, hr in raw_edges:
        pair = (old_to_new[ox_i], old_to_new[red_i])
        if pair in tree_pairs and id(hr) not in used_hr:
            bfs_hrs.append(hr)
            bfs_edges.append(pair)
            used_hr.add(id(hr))
            tree_pairs.discard(pair)
    for ox_i, red_i, hr in raw_edges:
        if id(hr) in used_hr:
            continue
        leftover_hrs.append(hr)
        leftover_edges.append((old_to_new[ox_i], old_to_new[red_i]))

    return oxidation_levels, bfs_hrs + leftover_hrs, bfs_edges + leftover_edges


def _extend_level_with_hydration(
    level: list[str],
    pka_pairs: Sequence[tuple[str, str, float]],
    hydration_edges: Sequence[HydrationEdge],
) -> list[str]:
    """Append hydrated forms and their acid-base chains to an oxidation level."""
    extended = list(level)
    changed = True
    while changed:
        changed = False
        for edge in hydration_edges:
            if edge.dehydrated not in extended or edge.hydrated in extended:
                continue
            extended.append(edge.hydrated)
            for name in _build_species_chain(edge.hydrated, pka_pairs)[1:]:
                if name not in extended:
                    extended.append(name)
            changed = True
    return extended


def build_pourbaix_graph(env) -> PourbaixGraph:
    """Build per-element redox chains and classified reaction edges from an environment."""
    half_reactions = sorted(
        getattr(env, "half_reactions", None) or [],
        key=lambda hr: hr.E0_SHE,
        reverse=True,
    )
    if not half_reactions:
        raise ValueError("Pourbaix graph requires at least one half-reaction.")

    acid_base, precipitation, association, hydration = _classify_all_reactions(env)
    pka_pairs = [(e.acid, e.base, e.pka) for e in acid_base]

    hr_by_element: dict[str, list] = {}
    for hr in half_reactions:
        ox = _primary_side_species(hr.oxidized, skip=POURBAIX_SKIP_SPECIES)
        red = _primary_side_species(hr.reduced, skip=POURBAIX_SKIP_SPECIES)
        el_ox = _redox_element(ox)
        el_red = _redox_element(red)
        if el_ox != el_red:
            raise ValueError(
                f"Half-reaction {hr.name or ox + '/' + red} spans elements {el_ox} and {el_red}; "
                "use separate chains per element."
            )
        hr_by_element.setdefault(el_ox, []).append(hr)

    chains: list[ElementChain] = []
    all_track: list[str] = []
    seen_track: set[str] = set()

    for element in sorted(hr_by_element):
        oxidation_levels, element_hrs, redox_edges = _build_oxidation_graph(
            element,
            hr_by_element[element],
            pka_pairs,
            hydration,
        )

        level_species = {name for level in oxidation_levels for name in level}
        chain_pka = [(a, b, p) for a, b, p in pka_pairs if a in level_species and b in level_species]
        chain_ab = [e for e in acid_base if e.acid in level_species and e.base in level_species]
        chain_assoc = [e for e in association if e.product in level_species or any(r in level_species for r in e.reactants)]
        chain_hyd = [e for e in hydration if e.dehydrated in level_species or e.hydrated in level_species]

        track: list[str] = []
        for level in oxidation_levels:
            for name in level:
                if name not in seen_track:
                    seen_track.add(name)
                    all_track.append(name)
                if name not in track:
                    track.append(name)
        for edge in chain_assoc:
            if edge.product not in track:
                track.append(edge.product)
                if edge.product not in seen_track:
                    seen_track.add(edge.product)
                    all_track.append(edge.product)
        for edge in precipitation:
            if edge.solid not in track:
                track.append(edge.solid)
                if edge.solid not in seen_track:
                    seen_track.add(edge.solid)
                    all_track.append(edge.solid)

        chains.append(
            ElementChain(
                element=element,
                oxidation_levels=oxidation_levels,
                half_reactions=element_hrs,
                track_species=track,
                pka_pairs=chain_pka,
                acid_base_edges=chain_ab,
                precipitation_edges=[e for e in precipitation],
                association_edges=chain_assoc,
                hydration_edges=chain_hyd,
                redox_edges=redox_edges,
            )
        )

    return PourbaixGraph(
        chains=chains,
        all_track_species=all_track,
        precipitation_edges=precipitation,
        acid_base_edges=acid_base,
        association_edges=association,
        hydration_edges=hydration,
    )


def _level_index(chain: ElementChain, species: str) -> int:
    for index, level in enumerate(chain.oxidation_levels):
        if species in level:
            return index
    raise ValueError(f"Species {species!r} not found in chain for {chain.element}.")


def _iter_redox_edges(chain: ElementChain) -> list[tuple[tuple[int, int], object]]:
    if chain.redox_edges and len(chain.redox_edges) == len(chain.half_reactions):
        return list(zip(chain.redox_edges, chain.half_reactions))
    return [((index, index + 1), hr) for index, hr in enumerate(chain.half_reactions)]


def _hr_connecting(chain: ElementChain, src: int, dst: int):
    for (ox_i, red_i), hr in _iter_redox_edges(chain):
        if ox_i == src and red_i == dst:
            return hr
    return None


def _best_path(chain: ElementChain, src: int, dst: int) -> list[tuple[int, int, object]]:
    """Shortest ox→red path, preferring atom-balanced half-reactions."""
    if src == dst:
        return []
    outgoing: dict[int, list[tuple[int, object]]] = defaultdict(list)
    for (ox_i, red_i), hr in _iter_redox_edges(chain):
        outgoing[ox_i].append((red_i, hr))

    best: dict[int, tuple[tuple[int, int, int], list[tuple[int, int, object]]]] = {
        src: ((0, 0, 0), [])
    }
    queue: deque[int] = deque([src])
    while queue:
        current = queue.popleft()
        cost, path = best[current]
        for nxt, hr in outgoing.get(current, []):
            unbalanced = 0 if _edge_atom_balanced(hr, chain.element) else 1
            new_path = path + [(current, nxt, hr)]
            new_cost = (cost[0] + unbalanced, cost[1] + 1, cost[2] + int(hr.n_electrons))
            prev = best.get(nxt)
            if prev is None or new_cost < prev[0]:
                best[nxt] = (new_cost, new_path)
                queue.append(nxt)
    found = best.get(dst)
    return found[1] if found is not None else []


def _prevalent_form(pH: float, level: Sequence[str], chain: ElementChain) -> str:
    dist = _within_level_distribution(pH, level, 1.0, chain)
    return max(level, key=lambda name: dist.get(name, 0.0))


def _hr_nominal_species(hr) -> tuple[str, str]:
    ox = _primary_side_species(hr.oxidized, skip=POURBAIX_SKIP_SPECIES)
    red = _primary_side_species(hr.reduced, skip=POURBAIX_SKIP_SPECIES)
    return ox, red


def _nominal_boundary_e(pH: float, hr, ox_chain: Sequence[str], red_chain: Sequence[str], pka_pairs, *, reduced_side: bool) -> float:
    if reduced_side:
        segments = _build_segments_reduced(hr, red_chain, pka_pairs)
    else:
        segments = _build_segments(hr, ox_chain, pka_pairs)
    return _boundary_e(segments, pH)


def _couple_Eh(
    pH: float,
    species_high: str,
    species_low: str,
    chain: ElementChain,
    hr,
) -> float:
    ox_chain = chain.oxidation_levels[_level_index(chain, species_high)]
    red_chain = chain.oxidation_levels[_level_index(chain, species_low)]
    nom_ox, nom_red = _hr_nominal_species(hr)
    n = hr.n_electrons
    m = hr.net_h_plus_stoichiometry()

    conc_ox = _within_level_distribution(pH, ox_chain, 1.0, chain)
    conc_red = _within_level_distribution(pH, red_chain, 1.0, chain)
    ox_high = max(conc_ox.get(species_high, 1e-300), 1e-300)
    ox_nom = max(conc_ox.get(nom_ox, 1e-300), 1e-300)
    red_low = max(conc_red.get(species_low, 1e-300), 1e-300)
    red_nom = max(conc_red.get(nom_red, 1e-300), 1e-300)

    ratio = (red_nom / red_low) / (ox_nom / ox_high)
    return hr.E0_SHE - (NERNST_K / n) * math.log10(max(ratio, 1e-300)) - (NERNST_K * m / n) * pH


def boundary_Eh(
    pH: float,
    species_high: str,
    species_low: str,
    chain: ElementChain,
) -> float:
    """
    Eh (V vs SHE) where ``species_high`` and ``species_low`` have equal formation energy at ``pH``.
    """
    pH = float(pH)
    level_high = _level_index(chain, species_high)
    level_low = _level_index(chain, species_low)
    if level_high == level_low:
        raise ValueError("species_high and species_low must be different oxidation levels.")
    if level_high > level_low:
        raise ValueError("species_high must be more oxidized than species_low.")

    path = _best_path(chain, level_high, level_low)
    direct = _hr_connecting(chain, level_high, level_low)
    if direct is not None:
        path = [(level_high, level_low, direct)]
    if not path:
        raise ValueError(
            f"No redox path from {species_high!r} to {species_low!r} in chain {chain.element}."
        )
    if len(path) == 1:
        return _couple_Eh(pH, species_high, species_low, chain, path[0][2])

    n_sum = 0.0
    nE_sum = 0.0
    for src, dst, hr in path:
        left = species_high if src == level_high else _prevalent_form(pH, chain.oxidation_levels[src], chain)
        right = species_low if dst == level_low else _prevalent_form(pH, chain.oxidation_levels[dst], chain)
        eh = _couple_Eh(pH, left, right, chain, hr)
        n_sum += hr.n_electrons
        nE_sum += hr.n_electrons * eh
    return nE_sum / n_sum


def _level_formation(chain: ElementChain, pH: float) -> tuple[list[str], list[float], list[float]]:
    """Prevalent forms plus (n_e per redox atom, E vs most-oxidized level) for each level."""
    forms = [_prevalent_form(pH, level, chain) for level in chain.oxidation_levels]
    n_levels = len(chain.oxidation_levels)
    n_e = [0.0] * n_levels
    e_vs_ref = [0.0] * n_levels
    ref = 0
    for index in range(1, n_levels):
        path = _best_path(chain, ref, index)
        direct = _hr_connecting(chain, ref, index)
        if direct is not None:
            path = [(ref, index, direct)]
        if not path:
            n_e[index] = float("inf")
            continue
        n_sum = 0.0
        nE_sum = 0.0
        for src, dst, hr in path:
            eh = _couple_Eh(pH, forms[src], forms[dst], chain, hr)
            n_sum += hr.n_electrons
            nE_sum += hr.n_electrons * eh
        n_atom = max(_count_element(chain.oxidation_levels[index][0], chain.element), 1)
        n_e[index] = n_sum / n_atom
        e_vs_ref[index] = nE_sum / n_sum
    return forms, n_e, e_vs_ref


def _select_hull_band(Eh: float, n_e: Sequence[float], e_vs_ref: Sequence[float]) -> int:
    """Pick the oxidation level with lowest formation energy per redox atom."""
    best = 0
    best_key = (0.0, 0.0, 0)
    for index, (n, e_ref) in enumerate(zip(n_e, e_vs_ref)):
        if n == float("inf"):
            continue
        energy = n * (Eh - e_ref)
        key = (energy, n, index)
        if index == 0 or key < best_key:
            best_key = key
            best = index
    return best


def _upper_hull_indices(n_e: Sequence[float], e_vs_ref: Sequence[float]) -> list[int]:
    """Level indices on the upper Frost hull of ``(n, n*E)``."""
    points: list[tuple[float, float, int]] = []
    for index, (n, e_ref) in enumerate(zip(n_e, e_vs_ref)):
        if n == float("inf"):
            continue
        points.append((n, n * e_ref, index))
    points.sort()
    merged: list[tuple[float, float, int]] = []
    for n, nE, index in points:
        if merged and abs(merged[-1][0] - n) < 1e-12:
            if nE > merged[-1][1] or (abs(nE - merged[-1][1]) < 1e-15 and index < merged[-1][2]):
                merged[-1] = (n, nE, index)
            continue
        merged.append((n, nE, index))

    def cross(
        origin: tuple[float, float, int],
        a: tuple[float, float, int],
        b: tuple[float, float, int],
    ) -> float:
        return (a[0] - origin[0]) * (b[1] - origin[1]) - (a[1] - origin[1]) * (b[0] - origin[0])

    hull: list[tuple[float, float, int]] = []
    for point in merged:
        while len(hull) >= 2 and cross(hull[-2], hull[-1], point) > 0:
            hull.pop()
        hull.append(point)
    return [index for _, _, index in hull]


def _hull_pair_eh(n_e: Sequence[float], e_vs_ref: Sequence[float], left: int, right: int) -> float:
    return (n_e[right] * e_vs_ref[right] - n_e[left] * e_vs_ref[left]) / (n_e[right] - n_e[left])


def _stable_hull_indices(
    n_e: Sequence[float],
    e_vs_ref: Sequence[float],
    *,
    min_width: float = 0.01,
) -> list[int]:
    """Upper hull with interior vertices removed when their Eh window is thinner than ``min_width``."""
    hull = _upper_hull_indices(n_e, e_vs_ref)
    while len(hull) > 2:
        bounds = [
            _hull_pair_eh(n_e, e_vs_ref, hull[index], hull[index + 1])
            for index in range(len(hull) - 1)
        ]
        drop: Optional[int] = None
        for index in range(1, len(hull) - 1):
            width = bounds[index - 1] - bounds[index]
            if width < min_width:
                drop = index
                break
        if drop is None:
            break
        del hull[drop]
    return hull


def _within_level_crossovers(
    level: Sequence[str],
    chain: ElementChain,
    *,
    pH_min: float = 0.0,
    pH_max: float = 14.0,
    steps: int = 400,
) -> list[tuple[str, str, float]]:
    """pH values where the prevalent species inside ``level`` switches."""
    if len(level) < 2:
        return []
    samples = np.linspace(pH_min, pH_max, max(steps, 2))
    pairs: list[tuple[str, str, float]] = []
    prev_name = _prevalent_form(float(samples[0]), level, chain)
    prev_pH = float(samples[0])
    for pH in samples[1:]:
        name = _prevalent_form(float(pH), level, chain)
        if name == prev_name:
            prev_pH = float(pH)
            continue
        lo, hi = prev_pH, float(pH)
        lo_name, hi_name = prev_name, name
        for _ in range(24):
            mid = 0.5 * (lo + hi)
            mid_name = _prevalent_form(mid, level, chain)
            if mid_name == lo_name:
                lo = mid
            else:
                hi = mid
                hi_name = mid_name
        pairs.append((lo_name, hi_name, 0.5 * (lo + hi)))
        prev_name = name
        prev_pH = float(pH)
    return pairs


def _solve_dimer(total: float, k: float) -> tuple[float, float]:
    """2A ⇌ A2 with C = [A] + 2[A2], K = [A2]/[A]^2."""
    if total <= 0 or k <= 0:
        return total, 0.0
    a = 2.0 * k
    b = 1.0
    c = -total
    disc = max(b * b - 4.0 * a * c, 0.0)
    a_conc = (-b + math.sqrt(disc)) / (2.0 * a)
    a_conc = max(a_conc, 0.0)
    dimer = k * a_conc * a_conc
    return a_conc, dimer


def _association_distribution(
    pH: float,
    level_species: Sequence[str],
    total: float,
    pka_pairs: Sequence[tuple[str, str, float]],
    association_edges: Sequence[AssociationEdge],
) -> dict[str, float]:
    base = _chain_concentrations(pH, level_species, pka_pairs, total)
    if total <= 0 or not association_edges:
        return base

    extra: dict[str, float] = {name: 0.0 for name in level_species}
    for edge in association_edges:
        reactant_names = list(edge.reactants)
        if not all(name in base or name in extra for name in reactant_names):
            continue
        if len(reactant_names) == 1 and edge.reactants[reactant_names[0]] == 2.0:
            monomer = reactant_names[0]
            monomer_total = base.get(monomer, 0.0)
            mono_c, dimer_c = _solve_dimer(monomer_total, edge.k)
            base[monomer] = mono_c
            extra[edge.product] = extra.get(edge.product, 0.0) + dimer_c
        elif len(reactant_names) == 2:
            a, b = reactant_names
            coeff_a = edge.reactants[a]
            coeff_b = edge.reactants[b]
            if coeff_a == 1.0 and coeff_b == 1.0:
                ca = max(base.get(a, 0.0), 1e-300)
                cb = max(base.get(b, 0.0), 1e-300)
                cab = edge.k * ca * cb
                base[a] = max(ca - cab, 1e-300)
                base[b] = max(cb - cab, 1e-300)
                extra[edge.product] = extra.get(edge.product, 0.0) + cab

    merged = dict(base)
    for name, value in extra.items():
        merged[name] = merged.get(name, 0.0) + value
    return merged


def _ksp_background(pH: float, extra: Optional[dict[str, float]] = None) -> dict[str, float]:
    background = {
        "H+": 10.0 ** (-float(pH)),
        "OH-": 10.0 ** (float(pH) - 14.0),
    }
    if extra:
        background.update(extra)
    return background


def _ksp_quotient(
    conc: dict[str, float],
    edge: PrecipitationEdge,
    background: Optional[dict[str, float]] = None,
) -> float:
    background = background or {}
    q = 1.0
    for ion, coeff in edge.ions.items():
        c = conc.get(ion, background.get(ion, 0.0))
        q *= max(c, 1e-300) ** coeff
    return q


def _check_ksp(
    conc: dict[str, float],
    precipitation_edges: Sequence[PrecipitationEdge],
    background: Optional[dict[str, float]] = None,
) -> Optional[str]:
    background = background or {}
    for edge in precipitation_edges:
        if _ksp_quotient(conc, edge, background) >= edge.ksp:
            return edge.solid
    return None


def _ksp_metal_ions(edge: PrecipitationEdge) -> list[str]:
    return [
        ion
        for ion in edge.ions
        if ion not in POURBAIX_SKIP_SPECIES and ion not in {"OH-", "H+"}
    ]


def _level_containing_ksp_ions(chain: ElementChain, edge: PrecipitationEdge) -> Optional[list[str]]:
    metals = _ksp_metal_ions(edge)
    if not metals:
        return None
    for level in chain.oxidation_levels:
        if any(ion in level for ion in metals):
            return list(level)
    return None


def _solid_from_ksp(
    pH: float,
    level: Sequence[str],
    c_tot: float,
    chain: ElementChain,
    precipitation_edges: Sequence[PrecipitationEdge],
    background: Optional[dict[str, float]] = None,
) -> Optional[str]:
    if c_tot <= 0:
        return None
    bg = _ksp_background(pH, background)
    conc = _within_level_distribution(pH, level, c_tot, chain)
    return _check_ksp(conc, precipitation_edges, bg)


def _ksp_crossover_pH(
    chain: ElementChain,
    edge: PrecipitationEdge,
    c_tot: float,
    *,
    pH_min: float,
    pH_max: float,
    steps: int = 400,
    background: Optional[dict[str, float]] = None,
) -> Optional[float]:
    level = _level_containing_ksp_ions(chain, edge)
    if level is None or c_tot <= 0:
        return None
    samples = np.linspace(pH_min, pH_max, max(int(steps), 2))
    prev_pH: Optional[float] = None
    prev_q: Optional[float] = None
    for pH in samples:
        pH_value = float(pH)
        conc = _within_level_distribution(pH_value, level, c_tot, chain)
        q = _ksp_quotient(conc, edge, _ksp_background(pH_value, background))
        if prev_q is not None and (prev_q - edge.ksp) * (q - edge.ksp) <= 0.0:
            if abs(q - prev_q) < 1e-300:
                return pH_value
            t = (edge.ksp - prev_q) / (q - prev_q)
            return float(prev_pH + t * (pH_value - prev_pH))
        prev_pH = pH_value
        prev_q = q
    return None


def _within_level_distribution(
    pH: float,
    level: Sequence[str],
    total: float,
    chain: ElementChain,
) -> dict[str, float]:
    for edge in chain.hydration_edges:
        if edge.dehydrated in level and edge.hydrated in level:
            sublevel = [edge.hydrated] + [
                name for name in level if name not in (edge.dehydrated, edge.hydrated)
            ]
            conc = _association_distribution(
                pH,
                sublevel,
                total,
                chain.pka_pairs,
                chain.association_edges,
            )
            hydrated_amount = conc.get(edge.hydrated, 0.0)
            conc[edge.dehydrated] = hydrated_amount / max(edge.kh, 1e-300)
            c_sum = sum(conc.get(name, 0.0) for name in level)
            if c_sum > 0:
                scale = total / c_sum
                for name in level:
                    conc[name] = conc.get(name, 0.0) * scale
            return conc

    return _association_distribution(
        pH,
        level,
        total,
        chain.pka_pairs,
        chain.association_edges,
    )


def _select_oxidation_band(Eh: float, bounds: Sequence[float]) -> int:
    """Backward-compatible window picker for monotonic adjacent bounds."""
    n_e = [0.0]
    e_vs_ref = [0.0]
    n_sum = 0.0
    nE_sum = 0.0
    for bound in bounds:
        n_sum += 1.0
        nE_sum += bound
        n_e.append(n_sum)
        e_vs_ref.append(nE_sum / n_sum)
    return _select_hull_band(Eh, n_e, e_vs_ref)


def dominant_at(
    pH: float,
    Eh: float,
    chain: ElementChain,
    c_tot: float,
    *,
    precipitation_edges: Optional[Sequence[PrecipitationEdge]] = None,
    background_ions: Optional[dict[str, float]] = None,
) -> dict[str, float]:
    levels = chain.oxidation_levels
    _, n_e, e_vs_ref = _level_formation(chain, pH)
    hull = _stable_hull_indices(n_e, e_vs_ref)
    hull_n = [n_e[index] if index in hull else float("inf") for index in range(len(n_e))]
    band = _select_hull_band(Eh, hull_n, e_vs_ref)

    conc = {name: 0.0 for name in chain.track_species}
    level_conc = _within_level_distribution(pH, levels[band], c_tot, chain)
    conc.update(level_conc)

    bg = dict(background_ions or {})
    bg.setdefault("H+", 10.0 ** (-float(pH)))
    bg.setdefault("OH-", 10.0 ** (float(pH) - 14.0))
    precip = precipitation_edges if precipitation_edges is not None else chain.precipitation_edges
    solid = _check_ksp(conc, precip, bg)
    if solid is not None:
        return {name: 0.0 for name in chain.track_species} | {solid: c_tot}

    return conc


def graph_speciation(
    pH: float,
    Eh: float,
    graph: PourbaixGraph,
    totals: dict[str, float],
    *,
    background_ions: Optional[dict[str, float]] = None,
) -> tuple[dict[str, float], int]:
    merged = {species: 0.0 for species in graph.all_track_species}
    for chain in graph.chains:
        total = totals.get(chain.element, 0.0)
        if total <= 0:
            continue
        chain_conc = dominant_at(
            pH,
            Eh,
            chain,
            total,
            precipitation_edges=graph.precipitation_edges,
            background_ions=background_ions,
        )
        for name, value in chain_conc.items():
            if name in merged:
                merged[name] = value
            elif name not in merged:
                merged[name] = value

    best_index = 0
    best_value = -1.0
    for index, name in enumerate(graph.all_track_species):
        value = merged.get(name, 0.0)
        if value > best_value:
            best_value = value
            best_index = index
    return merged, best_index


def boundary_species_at_pH(
    boundary: PourbaixBoundary,
    pH: float,
    graph: PourbaixGraph,
    totals: Optional[dict[str, float]] = None,
    background_ions: Optional[dict[str, float]] = None,
) -> tuple[str, str]:
    """Species in equilibrium along ``boundary`` at ``pH``."""
    if boundary.kind in {"acid_base", "ksp"}:
        return boundary.left, boundary.right
    if boundary.kind == "redox":
        for chain in graph.chains:
            if chain.element != boundary.chain_element:
                continue
            level_low, level_high = boundary.level_pair
            ox_level = chain.oxidation_levels[level_low]
            red_level = chain.oxidation_levels[level_high]
            left = _prevalent_form(pH, ox_level, chain)
            right = _prevalent_form(pH, red_level, chain)
            c_tot = float((totals or {}).get(chain.element, 0.0))
            solid = _solid_from_ksp(
                pH,
                ox_level,
                c_tot,
                chain,
                graph.precipitation_edges,
                background_ions,
            )
            if solid is not None:
                if left in ox_level:
                    left = solid
                if right in ox_level:
                    right = solid
            return left, right
    return boundary.left, boundary.right


def _interpolate_eh(boundary: PourbaixBoundary, pH: float) -> float:
    if len(boundary.pH) == 1:
        return float(boundary.Eh[0])
    return float(np.interp(pH, boundary.pH, boundary.Eh))


def _line_intersection(
    b1: PourbaixBoundary,
    b2: PourbaixBoundary,
) -> list[tuple[float, float]]:
    """Find intersections between two boundary polylines."""
    points: list[tuple[float, float]] = []
    for i in range(len(b1.pH) - 1):
        for j in range(len(b2.pH) - 1):
            x1, x2 = float(b1.pH[i]), float(b1.pH[i + 1])
            y1a, y2a = float(b1.Eh[i]), float(b1.Eh[i + 1])
            x3, x4 = float(b2.pH[j]), float(b2.pH[j + 1])
            y1b, y2b = float(b2.Eh[j]), float(b2.Eh[j + 1])

            denom = (x1 - x2) * (y1b - y2b) - (y1a - y2a) * (x3 - x4)
            if abs(denom) < 1e-15:
                continue
            t = ((x1 - x3) * (y1b - y2b) - (y1a - y1b) * (x3 - x4)) / denom
            u = -((x1 - x2) * (y1a - y1b) - (y1a - y2a) * (x1 - x3)) / denom
            if 0.0 <= t <= 1.0 and 0.0 <= u <= 1.0:
                px = x1 + t * (x2 - x1)
                py = y1a + t * (y2a - y1a)
                points.append((px, py))
    return points


def _cluster_points(points: list[tuple[float, float]], eps_pH: float, eps_eh: float) -> list[tuple[float, float]]:
    clusters: list[tuple[float, float]] = []
    for px, py in points:
        merged = False
        for index, (cx, cy) in enumerate(clusters):
            if abs(px - cx) <= eps_pH and abs(py - cy) <= eps_eh:
                clusters[index] = ((cx + px) / 2.0, (cy + py) / 2.0)
                merged = True
                break
        if not merged:
            clusters.append((px, py))
    return clusters


def compute_analytic_geometry(
    graph: PourbaixGraph,
    *,
    pH_min: float,
    pH_max: float,
    pH_steps: int = 200,
    eh_min: float = -1.5,
    eh_max: float = 1.5,
    junction_labels: Optional[dict[tuple[str, ...], str]] = None,
    totals: Optional[dict[str, float]] = None,
    background_ions: Optional[dict[str, float]] = None,
) -> tuple[list[PourbaixBoundary], list[PourbaixJunction]]:
    pH_samples = np.linspace(pH_min, pH_max, max(int(pH_steps), 2))
    boundaries: list[PourbaixBoundary] = []

    for chain in graph.chains:
        pair_points: dict[tuple[int, int], list[tuple[float, float, str, str]]] = defaultdict(list)
        pH_step = float(pH_samples[1] - pH_samples[0]) if len(pH_samples) > 1 else 0.05
        for pH in pH_samples:
            pH_value = float(pH)
            forms, n_e, e_vs_ref = _level_formation(chain, pH_value)
            hull = _stable_hull_indices(n_e, e_vs_ref)
            for left_i, right_i in zip(hull, hull[1:]):
                n_left = n_e[left_i]
                n_right = n_e[right_i]
                if abs(n_right - n_left) < 1e-15:
                    continue
                eh = (n_right * e_vs_ref[right_i] - n_left * e_vs_ref[left_i]) / (n_right - n_left)
                pair_points[(left_i, right_i)].append(
                    (pH_value, float(eh), forms[left_i], forms[right_i])
                )

        for (left_i, right_i), points in pair_points.items():
            runs: list[list[tuple[float, float, str, str]]] = []
            current: list[tuple[float, float, str, str]] = []
            for point in points:
                if current and point[0] - current[-1][0] > 1.5 * pH_step:
                    runs.append(current)
                    current = []
                current.append(point)
            if current:
                runs.append(current)
            for run in runs:
                if len(run) < 2:
                    continue
                boundaries.append(
                    PourbaixBoundary(
                        left=run[0][2],
                        right=run[0][3],
                        kind="redox",
                        pH=np.asarray([item[0] for item in run], dtype=float),
                        Eh=np.asarray([item[1] for item in run], dtype=float),
                        level_pair=(left_i, right_i),
                        chain_element=chain.element,
                    )
                )

        for level_index, level in enumerate(chain.oxidation_levels):
            for acid, base, pH_cross in _within_level_crossovers(level, chain, pH_min=pH_min, pH_max=pH_max):
                forms, n_e, e_vs_ref = _level_formation(chain, pH_cross)
                if level_index not in _stable_hull_indices(n_e, e_vs_ref):
                    continue
                boundaries.append(
                    PourbaixBoundary(
                        left=acid,
                        right=base,
                        kind="acid_base",
                        pH=np.asarray([pH_cross, pH_cross], dtype=float),
                        Eh=np.asarray([eh_min, eh_max], dtype=float),
                        level_pair=(level_index, level_index),
                        chain_element=chain.element,
                    )
                )

        if totals:
            c_tot = float(totals.get(chain.element, 0.0))
            for edge in graph.precipitation_edges:
                pH_cross = _ksp_crossover_pH(
                    chain,
                    edge,
                    c_tot,
                    pH_min=pH_min,
                    pH_max=pH_max,
                    steps=max(int(pH_steps) * 2, 400),
                    background=background_ions,
                )
                if pH_cross is None:
                    continue
                level = _level_containing_ksp_ions(chain, edge)
                if level is None:
                    continue
                aqueous = _prevalent_form(pH_cross, level, chain)
                boundaries.append(
                    PourbaixBoundary(
                        left=aqueous,
                        right=edge.solid,
                        kind="ksp",
                        pH=np.asarray([pH_cross, pH_cross], dtype=float),
                        Eh=np.asarray([eh_min, eh_max], dtype=float),
                        chain_element=chain.element,
                    )
                )

    boundaries.append(
        PourbaixBoundary(
            left="H2O",
            right="O2",
            kind="water",
            pH=pH_samples.copy(),
            Eh=1.229 - 2.0 * NERNST_K * pH_samples,
        )
    )
    boundaries.append(
        PourbaixBoundary(
            left="H+",
            right="H2",
            kind="water",
            pH=pH_samples.copy(),
            Eh=-2.0 * NERNST_K * pH_samples,
        )
    )

    raw_points: list[tuple[float, float, PourbaixBoundary, PourbaixBoundary]] = []
    for i, b1 in enumerate(boundaries):
        for b2 in boundaries[i + 1 :]:
            for px, py in _line_intersection(b1, b2):
                raw_points.append((px, py, b1, b2))

    point_map: dict[tuple[float, float], set[str]] = {}
    kind_map: dict[tuple[float, float], set[str]] = {}
    for px, py, b1, b2 in raw_points:
        key = (round(px, 4), round(py, 4))
        point_map.setdefault(key, set()).update({b1.left, b1.right, b2.left, b2.right})
        kind_map.setdefault(key, set()).update({b1.kind, b2.kind})

    clusters = _cluster_points([ (px, py) for px, py, _, _ in raw_points], eps_pH=0.05, eps_eh=0.02)
    junctions: list[PourbaixJunction] = []
    for cx, cy in clusters:
        species_set: set[str] = set()
        kinds: set[str] = set()
        for px, py, b1, b2 in raw_points:
            if abs(px - cx) <= 0.05 and abs(py - cy) <= 0.02:
                species_set.update({b1.left, b1.right, b2.left, b2.right})
                kinds.update({b1.kind, b2.kind})
        if len(species_set) < 3:
            continue
        species_tuple = tuple(sorted(species_set))
        label = format_junction_label(species_tuple)
        if junction_labels and species_tuple in junction_labels:
            label = junction_labels[species_tuple]
        junctions.append(
            PourbaixJunction(
                pH=cx,
                Eh=cy,
                species=species_tuple,
                label=label,
                boundary_kinds=tuple(sorted(kinds)),
            )
        )

    return boundaries, assign_junction_point_ids(junctions)


def element_totals_from_env(env, graph: PourbaixGraph) -> dict[str, float]:
    label_map = {compound.formula: float(env.concentrations[index]) for index, compound in enumerate(env.compounds)}
    totals: dict[str, float] = {}
    for chain in graph.chains:
        total = sum(label_map.get(name, 0.0) for name in chain.track_species)
        totals[chain.element] = total if total > 0 else 1.0
    return totals
