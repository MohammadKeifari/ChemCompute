"""Graph-based Pourbaix engine tests."""

from __future__ import annotations

import numpy as np
import pytest

from ChemCompute import Compound, Enviroment, HalfReaction, Pourbaix, Reaction
from ChemCompute._pourbaix_graph import (
    boundary_Eh,
    build_pourbaix_graph,
    format_junction_label,
    graph_speciation,
    element_totals_from_env,
)
from test_pourbaix import ka, selenium_environment


def aq(formula: str) -> Compound:
    from ChemCompute._formula import compound_from_species_token

    return compound_from_species_token(f"{formula}.aq")


def test_se_boundary_eh_spot_checks():
    chain = build_pourbaix_graph(selenium_environment()).chains[0]
    e_high = boundary_Eh(7.0, "SeO4-2", "SeO3-2", chain)
    e_low = boundary_Eh(7.0, "SeO3-2", "Se", chain)
    assert e_high == pytest.approx(0.50, abs=0.04)
    assert e_low == pytest.approx(0.26, abs=0.04)
    assert e_high > e_low


def test_se_seo3_dominance():
    env = selenium_environment()
    graph = build_pourbaix_graph(env)
    totals = element_totals_from_env(env, graph)

    _, idx_75 = graph_speciation(7.5, 0.35, graph, totals)
    assert graph.all_track_species[idx_75] == "SeO3-2"

    _, idx_10 = graph_speciation(10.0, 0.0, graph, totals)
    assert graph.all_track_species[idx_10] == "SeO3-2"


def test_fe_pourbaix_smoke():
    hr = HalfReaction.from_string(
        "Fe+3 & @e = Fe+2",
        concentrations=[0.01, 0.001],
        E0=0.771,
    )
    water = Compound("H2O", excess=True)
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
    env = Enviroment(
        kw,
        hr,
        concentrations={"Fe+3": 0.01, "Fe+2": 0.001},
        buffer=["H+"],
    )
    diagram = Pourbaix(
        env,
        pH_steps=5,
        Eh_steps=5,
        pH_min=1,
        pH_max=3,
        Eh_min=0.0,
        Eh_max=0.8,
    ).run()
    assert diagram.grid_dominant.shape == (5, 5)
    assert diagram.geometry_source == "analytic"
    assert len(diagram.analytic_boundaries) >= 1
    assert set(diagram.track_species) >= {"Fe+3", "Fe+2"}


def test_multi_element_chains():
    se_graph = build_pourbaix_graph(selenium_environment())
    assert len(se_graph.chains) == 1
    assert se_graph.chains[0].element == "Se"

    hr_i = HalfReaction.from_string(
        "I2 & @e = I-",
        concentrations=[1.0, 1.0],
        E0=0.54,
    )
    i_env = Enviroment(
        ka("HI", "I-", 10.0),
        hr_i,
        concentrations={"I2": 0.1},
        buffer=["H+"],
    )
    i_graph = build_pourbaix_graph(i_env)
    assert len(i_graph.chains) == 1
    assert i_graph.chains[0].element == "I"


def test_junction_label_format():
    label = format_junction_label(["SeO3-2", "HSeO3-", "SeO4-2"])
    assert label == "HSeO3- · SeO3-2 · SeO4-2"


def test_junction_point_ids_and_coords():
    diagram = Pourbaix(
        selenium_environment(),
        pH_min=0.0,
        pH_max=10.0,
        pH_steps=40,
        Eh_min=-1.0,
        Eh_max=1.2,
        Eh_steps=40,
    ).run()
    assert diagram.junction_points
    assert diagram.junction_points[0].point_id == "P1"
    assert diagram.analytic_junction_points
    table = diagram.junction_table(source="dominant")
    assert table[0]["point_id"] == "P1"
    pH, Eh = diagram.junction_coords("P1", source="dominant")
    assert pH == diagram.junction_points[0].pH


def test_boundary_modes():
    diagram = Pourbaix(
        selenium_environment(),
        pH_steps=25,
        Eh_steps=25,
        pH_min=0,
        pH_max=10,
    ).run()
    assert diagram.equal_boundary_lines
    assert diagram.analytic_boundaries
    assert len(diagram.analytic_junction_points) >= len(diagram.junction_points)


def test_junction_plot_label_styles():
    from ChemCompute._pourbaix_graph import PourbaixJunction, format_junction_plot_label

    junction = PourbaixJunction(
        pH=7.0,
        Eh=0.45,
        species=("SeO3-2", "SeO4-2", "HSeO3-"),
        label="HSeO3- · SeO3-2 · SeO4-2",
        point_id="P4",
    )
    assert format_junction_plot_label(junction, style="numbered") == "P4"
    assert format_junction_plot_label(junction, style="numbered_coords") == "P4(7.00, 0.45)"
    assert "SeO3-2" in format_junction_plot_label(junction, style="species")


def test_analytic_geometry_on_selenium():
    diagram = Pourbaix(
        selenium_environment(),
        pH_min=0.0,
        pH_max=10.0,
        pH_steps=20,
        Eh_min=-1.0,
        Eh_max=1.2,
        Eh_steps=20,
    ).run()
    assert diagram.geometry_source == "analytic"
    assert diagram.analytic_boundaries
    assert diagram.region_matrix().shape == (20, 20)
    assert diagram.dominant_species_at(7.5, 0.35) in diagram.track_species


def test_dominant_boundaries_include_acid_base_verticals():
    diagram = Pourbaix(
        selenium_environment(),
        pH_steps=40,
        Eh_steps=40,
        pH_min=0,
        pH_max=10,
        Eh_min=-1.0,
        Eh_max=1.2,
    ).run()
    track = diagram.track_species
    vertical_pairs = {
        frozenset({"H2SeO3", "HSeO3-"}),
        frozenset({"HSeO3-", "SeO3-2"}),
        frozenset({"H2Se", "HSe-"}),
    }
    drawn = {
        frozenset({track[left], track[right]})
        for left, right, pH_line, _ in diagram.equal_boundary_lines
        if len(pH_line) >= 2 and abs(float(pH_line[0]) - float(pH_line[-1])) < 1e-9
    }
    assert vertical_pairs <= drawn


def test_labeled_plot_style_smoke():
    import matplotlib

    matplotlib.use("Agg")
    diagram = Pourbaix(
        selenium_environment(),
        pH_steps=20,
        Eh_steps=20,
        pH_min=0,
        pH_max=10,
    ).run()
    ax = diagram.plot(plot_style="labeled", show=False)
    assert len(ax.texts) >= 3


def test_equilibrium_geometry_source():
    env = selenium_environment()
    diagram = Pourbaix(
        env,
        pH_steps=4,
        Eh_steps=4,
        speciation_method="equilibrium",
        max_iter=200,
    ).run()
    assert diagram.geometry_source == "grid"
    assert diagram.analytic_boundaries == []
