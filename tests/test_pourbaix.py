"""Pourbaix diagram tests and Selenium reference environment."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ChemCompute import Compound, Enviroment, HalfReaction, Pourbaix, Reaction
from ChemCompute._pourbaix_graph import build_pourbaix_graph, element_totals_from_env, graph_speciation

POURBAIX_OUTPUT_DIR = Path(__file__).resolve().parent.parent / "manual_test_output" / "pourbaix"


def aq(formula: str) -> Compound:
    from ChemCompute._formula import compound_from_species_token

    return compound_from_species_token(f"{formula}.aq")


def solid(formula: str) -> Compound:
    from ChemCompute._formula import compound_from_species_token

    return compound_from_species_token(f"{formula}.s")


def ka(acid: str, base: str, pka: float) -> Reaction:
    """Acid dissociation HA ⇌ H+ + A-. H+ concentration is managed by Pourbaix."""
    return Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": aq(acid), "rate_dependency": 1}],
        products=[
            {"stoichiometric_coefficient": 1, "compound": aq("H+"), "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": aq(base), "rate_dependency": 1},
        ],
        reactants_concentration=[0.0],
        products_concentration=[0.0, 0.0],
        K=10 ** (-pka),
    )


def selenium_environment(c_tot: float = 1.0) -> Enviroment:
    """Selenium Pourbaix reference system (C_tot Se = 1 M by default)."""
    water = Compound("H2O", excess=True, phase_point_list=[{"phase": "l", "temperature": 298}])
    kw = Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": water, "rate_dependency": 0}],
        products=[
            {"stoichiometric_coefficient": 1, "compound": aq("H+"), "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": aq("OH-"), "rate_dependency": 1},
        ],
        reactants_concentration=[0.0],
        products_concentration=[0.0, 0.0],
        K=1e-14,
    )

    hr1 = HalfReaction.from_string(
        "HSeO4- & 3_H+ & 2_@e = H2SeO3 & H2O.l",
        E0=1.15,
        name="HSeO4-/H2SeO3",
    )
    hr2 = HalfReaction.from_string(
        "H2SeO3 & 4_H+ & 4_@e = Se.s & 3_H2O.l",
        E0=0.74,
        name="H2SeO3/Se",
    )
    hr3 = HalfReaction.from_string(
        "Se.s & 2_H+ & 2_@e = H2Se",
        E0=-0.11,
        name="Se/H2Se",
    )

    return Enviroment(
        kw,
        ka("HSeO4-", "SeO4-2", 1.92),
        ka("H2SeO3", "HSeO3-", 2.62),
        ka("HSeO3-", "SeO3-2", 7.19),
        ka("H2Se", "HSe-", 3.89),
        hr1,
        hr2,
        hr3,
        concentrations={"HSeO4-": c_tot},
        buffer=["H+"],
    )


def carbon_environment(c_tot: float = 1e-3) -> Enviroment:
    """
    Carbon Pourbaix reference at 25 C with C_tot = 1 mM.

    Latimer diagram at pH 0 (V vs SHE):
    CO2 -0.52 V H2C2O4 -> HCO2H -0.03 V CH2O +0.13 V CH3OH
    with CO2 -0.20 V HCO2H branch (H2C2O4/HCO2H E0 = +0.32 V from Latimer).

    Acid-base: Kh(CO2/H2CO3)=1.7e-3; pKa H2CO3 6.4/10.3; H2C2O4 1.3/4.2; HCO2H 3.7.
    """
    water = Compound("H2O", excess=True, phase_point_list=[{"phase": "l", "temperature": 298}])
    kw = Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": water, "rate_dependency": 0}],
        products=[
            {"stoichiometric_coefficient": 1, "compound": aq("H+"), "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": aq("OH-"), "rate_dependency": 1},
        ],
        reactants_concentration=[0.0],
        products_concentration=[0.0, 0.0],
        K=1e-14,
    )
    hydration = Reaction(
        reactants=[
            {"stoichiometric_coefficient": 1, "compound": aq("CO2"), "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": water, "rate_dependency": 0},
        ],
        products=[{"stoichiometric_coefficient": 1, "compound": aq("H2CO3"), "rate_dependency": 1}],
        reactants_concentration=[0.0, 0.0],
        products_concentration=[0.0],
        K=1.7e-3,
    )

    hr_co2_oxalic = HalfReaction.from_string(
        "2_CO2 & 4_H+ & 2_@e = H2C2O4 & 2_H2O.l",
        E0=-0.52,
        name="CO2/H2C2O4",
    )
    hr_co2_formic = HalfReaction.from_string(
        "CO2 & 2_H+ & 2_@e = HCO2H & H2O.l",
        E0=-0.20,
        name="CO2/HCO2H",
    )
    hr_oxalic_formic = HalfReaction.from_string(
        "H2C2O4 & 2_H+ & 2_@e = HCO2H & H2O.l",
        E0=0.32,
        name="H2C2O4/HCO2H",
    )
    hr_formic_formaldehyde = HalfReaction.from_string(
        "HCO2H & 2_H+ & 2_@e = CH2O & H2O.l",
        E0=-0.03,
        name="HCO2H/CH2O",
    )
    hr_formaldehyde_methanol = HalfReaction.from_string(
        "CH2O & 2_H+ & 2_@e = CH3OH",
        E0=0.13,
        name="CH2O/CH3OH",
    )

    return Enviroment(
        kw,
        hydration,
        ka("H2CO3", "HCO3-", 6.4),
        ka("HCO3-", "CO3-2", 10.3),
        ka("H2C2O4", "HC2O4-", 1.3),
        ka("HC2O4-", "C2O4-2", 4.2),
        ka("HCO2H", "HCO2-", 3.7),
        hr_co2_oxalic,
        hr_co2_formic,
        hr_oxalic_formic,
        hr_formic_formaldehyde,
        hr_formaldehyde_methanol,
        concentrations={"CO2": c_tot},
        buffer=["H+"],
    )


def run_selenium_pourbaix(**kwargs) -> Pourbaix:
    defaults = {
        "pH_min": 0.0,
        "pH_max": 10.0,
        "pH_steps": 40,
        "Eh_min": -1.0,
        "Eh_max": 1.4,
        "Eh_steps": 40,
        "speciation_method": "model",
    }
    defaults.update(kwargs)
    return Pourbaix(selenium_environment(), **defaults).run()


def run_carbon_pourbaix(**kwargs) -> Pourbaix:
    defaults = {
        "pH_min": 0.0,
        "pH_max": 14.0,
        "pH_steps": 200,
        "Eh_min": -0.8,
        "Eh_max": 0.3,
        "Eh_steps": 200,
        "element_totals": {"C": 1e-3},
        "speciation_method": "model",
    }
    defaults.update(kwargs)
    return Pourbaix(carbon_environment(), **defaults).run()


def iron_environment(c_tot: float = 1e-3) -> Enviroment:
    """
    Iron Pourbaix at 25 C with Fe_tot = 1 mM.

    Latimer at pH 0: Fe+3 +0.77 V Fe+2 -0.44 V Fe(s).
    Hydrolysis pKa: Fe+3 2.2, Fe+2 9.5.
    Fe(OH)3(s) Ksp = [Fe+3][OH-]^3 = 1e-38
    (exam writes Fe+3 + 3 H2O ⇌ Fe(OH)3(s) + 3 H+ as that Ksp).
    """
    water = Compound("H2O", excess=True, phase_point_list=[{"phase": "l", "temperature": 298}])
    kw = Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": water, "rate_dependency": 0}],
        products=[
            {"stoichiometric_coefficient": 1, "compound": aq("H+"), "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": aq("OH-"), "rate_dependency": 1},
        ],
        reactants_concentration=[0.0],
        products_concentration=[0.0, 0.0],
        K=1e-14,
    )
    fe_oh3 = solid("Fe(OH)3")
    ksp = Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": fe_oh3, "rate_dependency": 1}],
        products=[
            {"stoichiometric_coefficient": 1, "compound": aq("Fe+3"), "rate_dependency": 1},
            {"stoichiometric_coefficient": 3, "compound": aq("OH-"), "rate_dependency": 1},
        ],
        reactants_concentration=[0.0],
        products_concentration=[0.0, 0.0],
        K=1e-38,
    )
    hr_fe3_fe2 = HalfReaction.from_string(
        "Fe+3 & @e = Fe+2",
        E0=0.77,
        name="Fe+3/Fe+2",
    )
    hr_fe2_fe = HalfReaction.from_string(
        "Fe+2 & 2_@e = Fe.s",
        E0=-0.44,
        name="Fe+2/Fe",
    )
    return Enviroment(
        kw,
        ksp,
        ka("Fe+3", "FeOH+2", 2.2),
        ka("Fe+2", "FeOH+", 9.5),
        hr_fe3_fe2,
        hr_fe2_fe,
        concentrations={"Fe+3": c_tot},
        buffer=["H+"],
    )


def run_iron_pourbaix(**kwargs) -> Pourbaix:
    defaults = {
        "pH_min": 0.0,
        "pH_max": 14.0,
        "pH_steps": 200,
        "Eh_min": -1.0,
        "Eh_max": 1.4,
        "Eh_steps": 200,
        "element_totals": {"Fe": 1e-3},
        "speciation_method": "model",
    }
    defaults.update(kwargs)
    return Pourbaix(iron_environment(), **defaults).run()


def amt_environment(c_tot: float = 1.0) -> Enviroment:
    """
    4-amino-TEMPO (AMT) Pourbaix at 25 C.

    Exam labels (not full piperidine structures):
      AMTO+ / AMTOH+2  oxoammonium O / amino-protonated O
      AMT / AMTH+      nitroxyl radical A / amino-protonated A
      AMTR / AMTRH+ / AMTRH2+2  hydroxylamine R and protonated forms

    Redox (reduction E°): AMTO+ / AMT = +0.74 V from A ⇌ O + e− (−0.74 V oxidation);
    AMT + H+ + e− ⇌ AMTR = +0.54 V as written.

    pKa = 14 − pKb: AMTH+ 9.0, AMTOH+2 6.0, AMTRH2+2 5.5, AMTRH+ 9.8.
    """
    water = Compound("H2O", excess=True, phase_point_list=[{"phase": "l", "temperature": 298}])
    kw = Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": water, "rate_dependency": 0}],
        products=[
            {"stoichiometric_coefficient": 1, "compound": aq("H+"), "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": aq("OH-"), "rate_dependency": 1},
        ],
        reactants_concentration=[0.0],
        products_concentration=[0.0, 0.0],
        K=1e-14,
    )
    hr_o_a = HalfReaction.from_string(
        "AMTO+ & @e = AMT",
        E0=0.74,
        name="AMTO+/AMT",
    )
    hr_a_r = HalfReaction.from_string(
        "AMT & H+ & @e = AMTR",
        E0=0.54,
        name="AMT/AMTR",
    )
    return Enviroment(
        kw,
        ka("AMTOH+2", "AMTO+", 6.0),
        ka("AMTH+", "AMT", 9.0),
        ka("AMTRH2+2", "AMTRH+", 5.5),
        ka("AMTRH+", "AMTR", 9.8),
        hr_o_a,
        hr_a_r,
        concentrations={"AMT": c_tot},
        buffer=["H+"],
    )


def run_amt_pourbaix(**kwargs) -> Pourbaix:
    defaults = {
        "pH_min": 0.0,
        "pH_max": 14.0,
        "pH_steps": 200,
        "Eh_min": -0.6,
        "Eh_max": 1.3,
        "Eh_steps": 200,
        "element_totals": {"A": 1.0},
        "speciation_method": "model",
    }
    defaults.update(kwargs)
    return Pourbaix(amt_environment(), **defaults).run()


@pytest.fixture(scope="module")
def se_diagram():
    return run_selenium_pourbaix()


def test_selenium_environment_track_species():
    graph = build_pourbaix_graph(selenium_environment())
    assert graph.chains[0].element == "Se"
    assert "SeO3-2" in graph.all_track_species
    assert "Se" in graph.all_track_species


def test_selenium_model_geometry(se_diagram):
    assert se_diagram.geometry_source == "analytic"
    assert se_diagram.analytic_boundaries
    assert se_diagram.junction_points
    assert se_diagram.equal_boundary_lines


def test_selenium_seo3_dominance_high_ph():
    env = selenium_environment()
    graph = build_pourbaix_graph(env)
    totals = element_totals_from_env(env, graph)
    _, idx = graph_speciation(7.5, 0.35, graph, totals)
    assert graph.all_track_species[idx] == "SeO3-2"


def test_selenium_vertical_acid_base_boundaries_drawn(se_diagram):
    track = se_diagram.track_species
    vertical_pairs = {
        frozenset({"H2SeO3", "HSeO3-"}),
        frozenset({"HSeO3-", "SeO3-2"}),
        frozenset({"H2Se", "HSe-"}),
    }
    drawn = {
        frozenset({track[left], track[right]})
        for left, right, pH_line, _ in se_diagram.equal_boundary_lines
        if len(pH_line) >= 2 and abs(float(pH_line[0]) - float(pH_line[-1])) < 1e-9
    }
    assert vertical_pairs <= drawn


def test_selenium_junction_table(se_diagram):
    table = se_diagram.junction_table(source="dominant")
    assert table[0]["point_id"] == "P1"
    pH, Eh = se_diagram.junction_coords("P1", source="dominant")
    assert pH == pytest.approx(se_diagram.junction_points[0].pH)


def test_selenium_labeled_plot_style():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    diagram = run_selenium_pourbaix(pH_steps=20, Eh_steps=20)
    ax = diagram.plot(plot_style="labeled", show=False)
    assert len(ax.texts) >= 3


def test_write_selenium_pourbaix_outputs():
    """Write reference Pourbaix figures under manual_test_output/pourbaix/."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    POURBAIX_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    diagram = run_selenium_pourbaix(pH_steps=200, Eh_steps=200)

    filled_path = POURBAIX_OUTPUT_DIR / "selenium_filled.png"
    labeled_path = POURBAIX_OUTPUT_DIR / "selenium_labeled.png"
    boundaries_path = POURBAIX_OUTPUT_DIR / "selenium_boundaries.png"

    diagram.plot(plot_style="filled", save=str(filled_path), show=False)
    diagram.plot(plot_style="labeled", save=str(labeled_path), show=False)
    diagram.plot_boundaries(save=str(boundaries_path), show=False)

    assert filled_path.is_file()
    assert labeled_path.is_file()
    assert boundaries_path.is_file()
    assert filled_path.stat().st_size > 0


def test_frame_intersections_on_selenium(se_diagram):
    points = se_diagram.frame_intersections()
    assert points
    table = se_diagram.frame_intersection_table()
    assert all("pH" in row and "Eh" in row and "edge" in row for row in table)
    edges = {row["edge"] for row in table}
    assert edges <= {"pH_min", "pH_max", "Eh_min", "Eh_max"}


def test_frame_intersections_plot_option():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    diagram = run_selenium_pourbaix(pH_steps=20, Eh_steps=20)
    ax = diagram.plot(show=False, show_frame_intersections=True, frame_intersection_labels=False)
    assert ax


def test_carbon_environment_track_species():
    graph = build_pourbaix_graph(carbon_environment())
    chain = graph.chains[0]
    assert chain.element == "C"
    assert chain.oxidation_levels[0][:4] == ["CO2", "H2CO3", "HCO3-", "CO3-2"]
    assert chain.oxidation_levels[1][:3] == ["H2C2O4", "HC2O4-", "C2O4-2"]
    formic_level = next(level for level in chain.oxidation_levels if "HCO2H" in level)
    assert formic_level == ["HCO2H", "HCO2-"]
    assert {hr.name for hr in chain.half_reactions} >= {
        "CO2/H2C2O4",
        "CO2/HCO2H",
        "H2C2O4/HCO2H",
        "HCO2H/CH2O",
        "CH2O/CH3OH",
    }
    assert ("CO2", "HCO2H") in {
        (
            chain.oxidation_levels[ox][0],
            chain.oxidation_levels[red][0],
        )
        for ox, red in chain.redox_edges
    }
    assert set(graph.all_track_species) >= {
        "CO2",
        "H2CO3",
        "HCO3-",
        "CO3-2",
        "H2C2O4",
        "HC2O4-",
        "C2O4-2",
        "HCO2H",
        "HCO2-",
        "CH2O",
        "CH3OH",
    }


def test_carbon_model_geometry():
    diagram = run_carbon_pourbaix(pH_steps=25, Eh_steps=25)
    assert diagram.geometry_source == "analytic"
    assert diagram.analytic_boundaries
    assert diagram.grid_dominant.shape == (25, 25)


def test_carbon_dominance_latimer_endpoints():
    diagram = run_carbon_pourbaix(pH_steps=30, Eh_steps=30)
    assert diagram.dominant_species_at(0.0, 0.2) == "CO2"
    assert diagram.dominant_species_at(0.0, -0.6) == "CH3OH"
    assert diagram.dominant_species_at(11.0, 0.0) in {"HCO3-", "CO3-2"}


def test_carbon_formic_region_mid_ph():
    graph = build_pourbaix_graph(carbon_environment())
    assert {"HCO2H", "H2C2O4", "HC2O4-", "CH2O"} <= set(graph.all_track_species)

    diagram = run_carbon_pourbaix(pH_steps=40, Eh_steps=40)
    from collections import Counter

    counts = Counter(
        diagram.track_species[diagram.grid_dominant[i, j]]
        for i in range(diagram.grid_dominant.shape[0])
        for j in range(diagram.grid_dominant.shape[1])
    )
    assert counts["CO2"] > 0
    assert counts["CH3OH"] > 0


def test_carbon_high_ph_carbonate_region():
    diagram = run_carbon_pourbaix(pH_steps=30, Eh_steps=30)
    assert diagram.dominant_species_at(11.0, 0.0) in {"HCO3-", "CO3-2"}


def test_carbon_carbonate_vertical_boundaries():
    diagram = run_carbon_pourbaix(pH_steps=40, Eh_steps=40)
    track = diagram.track_species
    verticals = []
    for left, right, pH_line, _ in diagram.equal_boundary_lines:
        if len(pH_line) >= 2 and abs(float(pH_line[0]) - float(pH_line[-1])) < 1e-9:
            verticals.append((frozenset({track[left], track[right]}), float(pH_line[0])))
    pairs = {pair for pair, _ in verticals}
    assert frozenset({"HCO3-", "CO3-2"}) in pairs
    assert frozenset({"CO2", "HCO3-"}) in pairs
    co2_hco3 = [pH for pair, pH in verticals if pair == frozenset({"CO2", "HCO3-"})]
    assert co2_hco3
    assert min(co2_hco3) == pytest.approx(9.17, abs=0.4)


def test_write_carbon_pourbaix_outputs():
    """Write carbon Pourbaix figures under manual_test_output/pourbaix/."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    POURBAIX_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    diagram = run_carbon_pourbaix(pH_steps=200, Eh_steps=200)

    filled_path = POURBAIX_OUTPUT_DIR / "carbon_filled.png"
    labeled_path = POURBAIX_OUTPUT_DIR / "carbon_labeled.png"
    boundaries_path = POURBAIX_OUTPUT_DIR / "carbon_boundaries.png"

    diagram.plot(plot_style="filled", save=str(filled_path), show=False)
    diagram.plot(plot_style="labeled", save=str(labeled_path), show=False)
    diagram.plot_boundaries(save=str(boundaries_path), show=False)

    assert filled_path.is_file()
    assert labeled_path.is_file()
    assert boundaries_path.is_file()
    assert filled_path.stat().st_size > 0


def test_iron_environment_track_species():
    graph = build_pourbaix_graph(iron_environment())
    chain = graph.chains[0]
    assert chain.element == "Fe"
    assert chain.oxidation_levels[0][:2] == ["Fe+3", "FeOH+2"]
    assert chain.oxidation_levels[1][:2] == ["Fe+2", "FeOH+"]
    assert "Fe" in chain.oxidation_levels[-1]
    assert {hr.name for hr in chain.half_reactions} >= {"Fe+3/Fe+2", "Fe+2/Fe"}
    assert "Fe(OH)3" in graph.all_track_species
    assert graph.precipitation_edges
    assert graph.precipitation_edges[0].ksp == pytest.approx(1e-38)
    assert "OH-" in graph.precipitation_edges[0].ions


def test_iron_dominance_latimer_and_hydroxide():
    diagram = run_iron_pourbaix(pH_steps=40, Eh_steps=40)
    assert diagram.dominant_species_at(0.0, 1.0) == "Fe+3"
    assert diagram.dominant_species_at(0.0, 0.2) == "Fe+2"
    assert diagram.dominant_species_at(0.0, -0.8) == "Fe"
    assert diagram.dominant_species_at(6.0, 1.0) == "Fe(OH)3"
    assert diagram.dominant_species_at(12.0, 0.8) == "Fe(OH)3"
    assert diagram.dominant_species_at(12.0, -0.2) in {"FeOH+", "Fe+2", "Fe"}


def test_iron_model_geometry():
    diagram = run_iron_pourbaix(pH_steps=25, Eh_steps=25)
    assert diagram.geometry_source == "analytic"
    assert diagram.analytic_boundaries
    assert diagram.grid_dominant.shape == (25, 25)


def test_iron_feoh3_boundary_lines():
    diagram = run_iron_pourbaix(pH_steps=40, Eh_steps=40)
    assert any(boundary.kind == "ksp" for boundary in diagram.analytic_boundaries)
    pairs = {
        frozenset({diagram.track_species[left], diagram.track_species[right]})
        for left, right, pH_line, eh_line in diagram.equal_boundary_lines
        if len(pH_line) >= 2
    }
    assert frozenset({"FeOH+2", "Fe(OH)3"}) in pairs
    assert frozenset({"Fe+2", "Fe(OH)3"}) in pairs
    assert frozenset({"FeOH+", "Fe(OH)3"}) in pairs
    p_alkaline = next(
        junction
        for junction in diagram.junction_points
        if {"Fe+2", "FeOH+", "Fe(OH)3"} <= set(junction.species)
    )
    assert p_alkaline.pH == pytest.approx(9.5, abs=0.08)


def test_write_iron_pourbaix_outputs():
    """Write iron Pourbaix figures under manual_test_output/pourbaix/."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    POURBAIX_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    diagram = run_iron_pourbaix(pH_steps=200, Eh_steps=200)
    species_sets = [set(junction.species) for junction in diagram.junction_points]
    assert any({"FeOH+2", "Fe+2", "Fe(OH)3"} <= names for names in species_sets)
    assert any({"Fe+3", "FeOH+2", "Fe+2"} <= names for names in species_sets)

    filled_path = POURBAIX_OUTPUT_DIR / "iron_filled.png"
    labeled_path = POURBAIX_OUTPUT_DIR / "iron_labeled.png"
    boundaries_path = POURBAIX_OUTPUT_DIR / "iron_boundaries.png"

    diagram.plot(plot_style="filled", save=str(filled_path), show=False)
    diagram.plot(plot_style="labeled", save=str(labeled_path), show=False)
    diagram.plot_boundaries(save=str(boundaries_path), show=False)

    assert filled_path.is_file()
    assert labeled_path.is_file()
    assert boundaries_path.is_file()
    assert filled_path.stat().st_size > 0


def test_amt_environment_track_species():
    graph = build_pourbaix_graph(amt_environment())
    chain = graph.chains[0]
    assert chain.element == "A"
    assert chain.oxidation_levels[0] == ["AMTOH+2", "AMTO+"]
    assert chain.oxidation_levels[1] == ["AMTH+", "AMT"]
    assert chain.oxidation_levels[2] == ["AMTRH2+2", "AMTRH+", "AMTR"]
    assert {hr.name for hr in chain.half_reactions} >= {"AMTO+/AMT", "AMT/AMTR"}
    assert set(graph.all_track_species) >= {
        "AMTO+",
        "AMTOH+2",
        "AMT",
        "AMTH+",
        "AMTR",
        "AMTRH+",
        "AMTRH2+2",
    }


def test_amt_dominance_latimer_and_pka():
    diagram = run_amt_pourbaix(pH_steps=40, Eh_steps=40)
    assert diagram.dominant_species_at(0.0, 1.2) == "AMTOH+2"
    # At pH 0 the AMTH+ band is only ~5 mV (O/A ≈ A/R ≈ 0.91 V).
    assert diagram.dominant_species_at(0.0, 0.8) == "AMTRH2+2"
    assert diagram.dominant_species_at(0.0, 0.0) == "AMTRH2+2"
    assert diagram.dominant_species_at(2.0, 0.8) == "AMTH+"
    assert diagram.dominant_species_at(7.0, 1.1) == "AMTO+"
    assert diagram.dominant_species_at(7.0, 0.55) == "AMTH+"
    assert diagram.dominant_species_at(11.0, 1.0) == "AMTO+"
    assert diagram.dominant_species_at(11.0, 0.5) == "AMT"
    assert diagram.dominant_species_at(11.0, -0.3) == "AMTR"


def test_amt_model_geometry():
    diagram = run_amt_pourbaix(pH_steps=25, Eh_steps=25)
    assert diagram.geometry_source == "analytic"
    assert diagram.analytic_boundaries
    assert diagram.junction_points
    assert diagram.grid_dominant.shape == (25, 25)


def _amt_junction(diagram, names: set[str]) -> dict:
    target = frozenset(names)
    for row in diagram.junction_table():
        if frozenset(row["species"]) == target:
            return row
    raise AssertionError(f"No dominant junction for {sorted(names)}")


def test_amt_triple_points():
    diagram = run_amt_pourbaix(pH_steps=50, Eh_steps=50)
    p_r1 = _amt_junction(diagram, {"AMTH+", "AMTRH+", "AMTRH2+2"})
    p_o = _amt_junction(diagram, {"AMTH+", "AMTO+", "AMTOH+2"})
    p_a_redox = _amt_junction(diagram, {"AMT", "AMTH+", "AMTO+"})
    p_a_ar = _amt_junction(diagram, {"AMT", "AMTH+", "AMTRH+"})
    p_r2 = _amt_junction(diagram, {"AMT", "AMTR", "AMTRH+"})
    assert p_r1["pH"] == pytest.approx(5.5, abs=0.3)
    assert p_o["pH"] == pytest.approx(6.0, abs=0.3)
    assert p_o["Eh"] == pytest.approx(0.74 + 0.05916 * 3.0, abs=0.08)
    assert p_a_redox["pH"] == pytest.approx(9.0, abs=0.3)
    assert p_a_redox["Eh"] == pytest.approx(0.74, abs=0.08)
    assert p_a_ar["pH"] == pytest.approx(9.0, abs=0.3)
    assert p_r2["pH"] == pytest.approx(9.8, abs=0.3)


def test_write_amt_pourbaix_outputs():
    """Write AMT Pourbaix figures under manual_test_output/pourbaix/."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    POURBAIX_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    diagram = run_amt_pourbaix(pH_steps=200, Eh_steps=200)

    filled_path = POURBAIX_OUTPUT_DIR / "amt_filled.png"
    labeled_path = POURBAIX_OUTPUT_DIR / "amt_labeled.png"
    boundaries_path = POURBAIX_OUTPUT_DIR / "amt_boundaries.png"

    diagram.plot(plot_style="filled", save=str(filled_path), show=False)
    diagram.plot(plot_style="labeled", save=str(labeled_path), show=False)
    diagram.plot_boundaries(save=str(boundaries_path), show=False)

    assert filled_path.is_file()
    assert labeled_path.is_file()
    assert boundaries_path.is_file()
    assert filled_path.stat().st_size > 0


if __name__ == "__main__":
    test_write_selenium_pourbaix_outputs()
    test_write_carbon_pourbaix_outputs()
    test_write_iron_pourbaix_outputs()
    test_write_amt_pourbaix_outputs()
    print(f"Wrote Pourbaix plots to {POURBAIX_OUTPUT_DIR.resolve()}")
