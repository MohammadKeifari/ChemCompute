"""Formula charge inference and unified reaction string notation."""

from __future__ import annotations

import pytest

from ChemCompute import Enviroment, HalfReaction, Reaction
from ChemCompute._formula import (
    compound_from_species_token,
    infer_ionic_charge,
    normalize_formula,
)


@pytest.mark.parametrize(
    "formula,expected",
    [
        ("H+", 1),
        ("OH-", -1),
        ("NH4+", 1),
        ("Fe+3", 3),
        ("SeO4-2", -2),
        ("Fe(CN)6-4", -4),
        ("[Fe(CN)6]-4", -4),
        ("H2O", 0),
        ("Se.s", 0),
    ],
)
def test_infer_ionic_charge(formula, expected):
    assert infer_ionic_charge(formula) == expected


@pytest.mark.parametrize(
    "name,canonical",
    [
        ("[Fe(CN)6]-4", "Fe(CN)6-4"),
        ("Fe(CN)6-4", "Fe(CN)6-4"),
        ("H+.aq", "H+.aq"),
        ("[Fe(CN)6]-4.aq", "Fe(CN)6-4.aq"),
    ],
)
def test_normalize_formula(name, canonical):
    assert normalize_formula(name) == canonical


def test_compound_from_species_token_sets_charge_and_phase():
    compound = compound_from_species_token("SeO4-2.aq")
    assert compound.formula == "SeO4-2"
    assert compound.charge == -2
    assert compound.phase_point_list == [{"temperature": 298, "phase": "aq"}]


def test_reaction_from_string_sets_compound_charge():
    rxn = Reaction.from_string("H2SeO3 > H+ & HSeO3-", concentrations=[0.0, 0.0, 0.0], K=10 ** -2.62)
    charges = {entry["compound"].formula: entry["compound"].charge for entry in rxn.compounds}
    assert charges["H+"] == 1
    assert charges["HSeO3-"] == -1
    assert charges["H2SeO3"] == 0


def test_half_reaction_from_string_sets_compound_charge():
    hr = HalfReaction.from_string(
        "HSeO4- & 3_H+ & 2_@e = H2SeO3 & H2O.l",
        concentrations=[1.0, 0.0, 0.0, 0.0],
        E0=1.15,
    )
    formulas = {entry["compound"].formula for entry in hr.compounds}
    assert "HSeO4-" in formulas
    assert any(entry["compound"].formula == "HSeO4-" and entry["compound"].charge == -1 for entry in hr.compounds)


def test_environment_auto_fills_charge_map():
    kw = Reaction.from_string("H2O.l > H+ & OH-", K=1e-14, concentrations=[0.0, 0.0, 0.0])
    ka = Reaction.from_string("H2SeO3 > H+ & HSeO3-", K=10 ** -2.62, concentrations=[0.0, 0.0, 0.0])
    hr = HalfReaction.from_string(
        "HSeO4- & 3_H+ & 2_@e = H2SeO3 & H2O.l",
        concentrations=[1.0, 0.0, 0.0, 0.0],
        E0=1.15,
    )
    env = Enviroment(kw, ka, hr, concentrations={"HSeO4-": 1.0}, buffer=["H+"])
    assert env.charge_map["H+"] == 1
    assert env.charge_map["OH-"] == -1
    assert env.charge_map["HSeO3-"] == -1
    assert env.charge_map["HSeO4-"] == -1
