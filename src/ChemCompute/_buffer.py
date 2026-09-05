"""Buffer capacity and Henderson-Hasselbalch diagnostics."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from ._activity import resolve_charge


@dataclass
class BufferPairDiagnostics:
    """Diagnostics for one conjugate acid/base pair."""

    acid_formula: str
    base_formula: str
    pka: float
    acid_concentration: float
    base_concentration: float
    pair_beta: float
    hh_predicted_pH: float
    hh_error: float
    hh_valid: bool


@dataclass
class BufferDiagnostics:
    """Buffer capacity and Henderson-Hasselbalch summary."""

    pH: float
    beta: float
    buffer_pairs: list[BufferPairDiagnostics] = field(default_factory=list)
    hh_predictions: dict[str, float] = field(default_factory=dict)


def _species_activity(env, formula: str, concentration: float, concentrations: np.ndarray) -> float:
    if getattr(env, "activity_model", None) is None:
        return concentration
    for j, compound in enumerate(env.compounds):
        if compound.formula == formula:
            gammas = env.activity_model.gamma_array(env, concentrations, env.T)
            return gammas[j] * concentration
    return concentration


def _find_h_plus_index(env) -> Optional[int]:
    for j, compound in enumerate(env.compounds):
        if compound.formula in ("H+", "H+"):
            return j
    for j, compound in enumerate(env.compounds):
        if compound.formula.replace(" ", "") in ("H+", "H+"):
            return j
    return None


def _identify_buffer_pairs(env):
    """Find HA + H+ -> A- style equilibria (single proton transfer)."""
    pairs = []
    h_index = _find_h_plus_index(env)
    if h_index is None:
        return pairs
    h_formula = env.compounds[h_index].formula

    for rxn in env.reactions:
        if getattr(rxn, "infinite_K", False):
            continue
        reactants = rxn.reactants
        products = rxn.products
        reactant_formulas = {r["compound"].formula for r in reactants}
        product_formulas = {p["compound"].formula for p in products}
        if h_formula not in reactant_formulas and h_formula not in product_formulas:
            continue

        if h_formula in reactant_formulas and len(reactants) == 2 and len(products) == 1:
            other_r = [r for r in reactants if r["compound"].formula != h_formula][0]
            acid = products[0]["compound"].formula
            base = other_r["compound"].formula
        elif h_formula in product_formulas and len(reactants) == 1 and len(products) == 2:
            other_p = [p for p in products if p["compound"].formula != h_formula][0]
            acid = reactants[0]["compound"].formula
            base = other_p["compound"].formula
        elif len(reactants) == 2 and len(products) == 2:
            if h_formula in reactant_formulas:
                other_r = [r for r in reactants if r["compound"].formula != h_formula][0]
                other_p = [p for p in products if p["compound"].formula != h_formula][0]
                acid = other_p["compound"].formula
                base = other_r["compound"].formula
            else:
                other_r = [r for r in reactants if r["compound"].formula != h_formula][0]
                other_p = [p for p in products if p["compound"].formula != h_formula][0]
                acid = other_r["compound"].formula
                base = other_p["compound"].formula
        else:
            continue

        ka = rxn.K
        pka = -math.log10(max(ka, 1e-300))
        pairs.append(
            {
                "acid": acid,
                "base": base,
                "ka": ka,
                "pka": pka,
            }
        )
    return pairs


def buffer_diagnostics(env, equilibrium_concentrations=None) -> BufferDiagnostics:
    """
    Compute buffer capacity beta(pH) and Henderson-Hasselbalch diagnostics.

    beta = 2.303 * sum(Ka * [HA]) over identified weak-acid pairs.
    """
    if equilibrium_concentrations is None:
        concentrations = np.array(env.concentrations, dtype=float)
    else:
        concentrations = np.array(equilibrium_concentrations, dtype=float)

    h_index = _find_h_plus_index(env)
    if h_index is None:
        return BufferDiagnostics(pH=float("nan"), beta=0.0)

    h_conc = concentrations[h_index]
    h_activity = _species_activity(env, env.compounds[h_index].formula, h_conc, concentrations)
    pH = -math.log10(max(h_activity, 1e-300))

    conc_dict = {env.compounds[j].formula: concentrations[j] for j in range(len(env.compounds))}
    pairs_info = _identify_buffer_pairs(env)
    buffer_pairs = []
    beta = 0.0
    hh_predictions = {}

    for pair in pairs_info:
        acid_c = conc_dict.get(pair["acid"], 0.0)
        base_c = conc_dict.get(pair["base"], 0.0)
        acid_a = _species_activity(env, pair["acid"], acid_c, concentrations)
        base_a = _species_activity(env, pair["base"], base_c, concentrations)
        ka = pair["ka"]
        pair_beta = 2.303 * ka * acid_a
        beta += pair_beta

        if acid_a > 0 and base_a > 0:
            hh_pH = pair["pka"] + math.log10(base_a / acid_a)
        else:
            hh_pH = float("nan")
        hh_error = abs(hh_pH - pH) if not math.isnan(hh_pH) else float("nan")
        hh_valid = hh_error < 0.5 if not math.isnan(hh_error) else False

        buffer_pairs.append(
            BufferPairDiagnostics(
                acid_formula=pair["acid"],
                base_formula=pair["base"],
                pka=pair["pka"],
                acid_concentration=acid_c,
                base_concentration=base_c,
                pair_beta=pair_beta,
                hh_predicted_pH=hh_pH,
                hh_error=hh_error,
                hh_valid=hh_valid,
            )
        )
        hh_predictions[pair["acid"]] = hh_pH

    return BufferDiagnostics(
        pH=pH,
        beta=beta,
        buffer_pairs=buffer_pairs,
        hh_predictions=hh_predictions,
    )
