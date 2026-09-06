"""Graph-based Pourbaix speciation, boundaries, and analytic geometry."""

from __future__ import annotations

import math
import re
import warnings
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
    chain = [start]
    current = start
    pair_map = {acid: base for acid, base, _ in pka_pairs}
    while current in pair_map:
        nxt = pair_map[current]
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
        pka = pair_map[(species[index], species[index + 1])]
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
            raise ValueError(f"Missing pKa pair for {acid} -> {base}.")
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
            raise ValueError(f"Missing pKa pair for {acid} -> {base}.")
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
    aq_r = [e for e in reactants if not _is_solid(e) and _formula(e) not in POURBAIX_SKIP_SPECIES]
    aq_p = [e for e in products if not _is_solid(e) and _formula(e) not in POURBAIX_SKIP_SPECIES]

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
        element_hrs = sorted(hr_by_element[element], key=lambda hr: hr.E0_SHE, reverse=True)
        oxidation_levels: list[list[str]] = []

        first_ox = _primary_side_species(element_hrs[0].oxidized, skip=POURBAIX_SKIP_SPECIES)
        oxidation_levels.append(_build_species_chain(first_ox, pka_pairs))

        for index in range(len(element_hrs) - 1):
            bridge = _primary_side_species(element_hrs[index].reduced, skip=POURBAIX_SKIP_SPECIES)
            oxidation_levels.append(_build_species_chain(bridge, pka_pairs))

        last_red = _primary_side_species(element_hrs[-1].reduced, skip=POURBAIX_SKIP_SPECIES)
        oxidation_levels.append(_build_species_chain(last_red, pka_pairs))

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


def _prevalent_form(pH: float, level: Sequence[str], pka_pairs: Sequence[tuple[str, str, float]]) -> str:
    return _acid_base_form(pH, level, pka_pairs)


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


def boundary_Eh(
    pH: float,
    species_high: str,
    species_low: str,
    chain: ElementChain,
) -> float:
    """
    Eh (V vs SHE) where ``species_high`` and ``species_low`` have equal concentration at ``pH``.
    """
    pH = float(pH)
    level_high = _level_index(chain, species_high)
    level_low = _level_index(chain, species_low)
    if level_high >= level_low:
        raise ValueError("species_high must be more oxidized than species_low.")

    if level_low - level_high != 1:
        mid = _prevalent_form(pH, chain.oxidation_levels[level_high + 1], chain.pka_pairs)
        return boundary_Eh(pH, species_high, mid, chain)

    hr = chain.half_reactions[level_high]
    ox_chain = chain.oxidation_levels[level_high]
    red_chain = chain.oxidation_levels[level_low]
    nom_ox, nom_red = _hr_nominal_species(hr)
    n = hr.n_electrons
    m = hr.net_h_plus_stoichiometry()

    conc_ox = _chain_concentrations(pH, ox_chain, chain.pka_pairs, 1.0)
    conc_red = _chain_concentrations(pH, red_chain, chain.pka_pairs, 1.0)
    ox_high = max(conc_ox.get(species_high, 1e-300), 1e-300)
    ox_nom = max(conc_ox.get(nom_ox, 1e-300), 1e-300)
    red_low = max(conc_red.get(species_low, 1e-300), 1e-300)
    red_nom = max(conc_red.get(nom_red, 1e-300), 1e-300)

    ratio = (red_nom / red_low) / (ox_nom / ox_high)
    return hr.E0_SHE - (NERNST_K / n) * math.log10(max(ratio, 1e-300)) - (NERNST_K * m / n) * pH


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


def _check_ksp(
    conc: dict[str, float],
    precipitation_edges: Sequence[PrecipitationEdge],
    background: Optional[dict[str, float]] = None,
) -> Optional[str]:
    background = background or {}
    for edge in precipitation_edges:
        q = 1.0
        for ion, coeff in edge.ions.items():
            c = conc.get(ion, background.get(ion, 0.0))
            q *= max(c, 1e-300) ** coeff
        if q >= edge.ksp:
            return edge.solid
    return None


def _within_level_distribution(
    pH: float,
    level: Sequence[str],
    total: float,
    chain: ElementChain,
) -> dict[str, float]:
    conc = _association_distribution(
        pH,
        level,
        total,
        chain.pka_pairs,
        chain.association_edges,
    )
    for edge in chain.hydration_edges:
        if edge.dehydrated in level and edge.hydrated in level:
            cd = conc.get(edge.dehydrated, 0.0)
            ch = conc.get(edge.hydrated, 0.0)
            if cd + ch > 0:
                kh = edge.kh
                ch_new = kh * cd
                scale = total / max(cd + ch_new, 1e-300)
                conc[edge.dehydrated] = cd * scale
                conc[edge.hydrated] = ch_new * scale
    return conc


def _near_boundary(Eh: float, bounds: Sequence[float], band: int, width: float = 0.015) -> Optional[int]:
    if band > 0 and -width < Eh - bounds[band - 1] < width:
        return band - 1
    if band < len(bounds) and -width < Eh - bounds[band] < width:
        return band + 1
    return None


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
    n = len(levels)
    forms = [_prevalent_form(pH, level, chain.pka_pairs) for level in levels]
    bounds = [boundary_Eh(pH, forms[i], forms[i + 1], chain) for i in range(n - 1)]

    band = n - 1
    for index, bound in enumerate(bounds):
        if Eh >= bound:
            band = index
            break

    conc = {name: 0.0 for name in chain.track_species}
    level_conc = _within_level_distribution(pH, levels[band], c_tot, chain)
    conc.update(level_conc)

    neighbor = _near_boundary(Eh, bounds, band)
    if neighbor is not None and 0 <= neighbor < n:
        other = _within_level_distribution(pH, levels[neighbor], c_tot * 0.5, chain)
        for name in chain.track_species:
            conc[name] = level_conc.get(name, 0.0) * 0.5 + other.get(name, 0.0) * 0.5

    precip = precipitation_edges if precipitation_edges is not None else chain.precipitation_edges
    solid = _check_ksp(conc, precip, background_ions)
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


def boundary_species_at_pH(boundary: PourbaixBoundary, pH: float, graph: PourbaixGraph) -> tuple[str, str]:
    """Species in equilibrium along ``boundary`` at ``pH``."""
    if boundary.kind == "acid_base":
        return boundary.left, boundary.right
    if boundary.kind == "redox":
        for chain in graph.chains:
            if chain.element != boundary.chain_element:
                continue
            level_low, level_high = boundary.level_pair
            left = _prevalent_form(pH, chain.oxidation_levels[level_low], chain.pka_pairs)
            right = _prevalent_form(pH, chain.oxidation_levels[level_high], chain.pka_pairs)
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
) -> tuple[list[PourbaixBoundary], list[PourbaixJunction]]:
    pH_samples = np.linspace(pH_min, pH_max, max(int(pH_steps), 2))
    boundaries: list[PourbaixBoundary] = []

    for chain in graph.chains:
        n = len(chain.oxidation_levels)
        for i in range(n - 1):
            eh_curve = []
            left_labels = []
            right_labels = []
            for pH in pH_samples:
                left = _prevalent_form(float(pH), chain.oxidation_levels[i], chain.pka_pairs)
                right = _prevalent_form(float(pH), chain.oxidation_levels[i + 1], chain.pka_pairs)
                left_labels.append(left)
                right_labels.append(right)
                eh_curve.append(boundary_Eh(float(pH), left, right, chain))
            boundaries.append(
                PourbaixBoundary(
                    left=left_labels[0],
                    right=right_labels[0],
                    kind="redox",
                    pH=pH_samples.copy(),
                    Eh=np.asarray(eh_curve, dtype=float),
                    level_pair=(i, i + 1),
                    chain_element=chain.element,
                )
            )

        seen_pka: set[tuple[str, str]] = set()
        for level_index, level in enumerate(chain.oxidation_levels):
            for acid, base, pka in chain.pka_pairs:
                if (acid, base) in seen_pka:
                    continue
                if acid in level and base in level:
                    seen_pka.add((acid, base))
                    boundaries.append(
                        PourbaixBoundary(
                            left=acid,
                            right=base,
                            kind="acid_base",
                            pH=np.asarray([pka, pka], dtype=float),
                            Eh=np.asarray([eh_min, eh_max], dtype=float),
                            level_pair=(level_index, level_index),
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
