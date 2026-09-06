"""Pourbaix diagram tests and Selenium reference environment."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ChemCompute import Compound, Enviroment, HalfReaction, Pourbaix, Reaction
from ChemCompute._pourbaix_graph import build_pourbaix_graph, element_totals_from_env, graph_speciation

POURBAIX_OUTPUT_DIR = Path(__file__).resolve().parent.parent / "manual_test_output" / "pourbaix"


def aq(formula: str, *, charge: int = 0) -> Compound:
    return Compound(
        formula,
        phase_point_list=[{"phase": "aq", "temperature": 298}],
        charge=charge,
    )


def ka(acid: str, base: str, pka: float, *, charge_acid: int = 0, charge_base: int = -1) -> Reaction:
    """Acid dissociation HA ⇌ H+ + A-. H+ concentration is managed by Pourbaix."""
    return Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": aq(acid, charge=charge_acid), "rate_dependency": 1}],
        products=[
            {"stoichiometric_coefficient": 1, "compound": aq("H+", charge=1), "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": aq(base, charge=charge_base), "rate_dependency": 1},
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
            {"stoichiometric_coefficient": 1, "compound": aq("H+", charge=1), "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": aq("OH-", charge=-1), "rate_dependency": 1},
        ],
        reactants_concentration=[0.0],
        products_concentration=[0.0, 0.0],
        K=1e-14,
    )

    hr1 = HalfReaction.from_string_complex_syntax(
        "HSeO4- & 3_H+ + 2_@e = H2SeO3 & H2O.l",
        E0=1.15,
        name="HSeO4-/H2SeO3",
    )
    hr2 = HalfReaction.from_string_complex_syntax(
        "H2SeO3 & 4_H+ + 4_@e = Se.s & 3_H2O.l",
        E0=0.74,
        name="H2SeO3/Se",
    )
    hr3 = HalfReaction.from_string_complex_syntax(
        "Se.s & 2_H+ + 2_@e = H2Se",
        E0=-0.11,
        name="Se/H2Se",
    )

    return Enviroment(
        kw,
        ka("HSeO4-", "SeO4-2", 1.92, charge_acid=-1, charge_base=-2),
        ka("H2SeO3", "HSeO3-", 2.62),
        ka("HSeO3-", "SeO3-2", 7.19, charge_acid=-1, charge_base=-2),
        ka("H2Se", "HSe-", 3.89),
        hr1,
        hr2,
        hr3,
        concentrations={"HSeO4-": c_tot},
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
    diagram = run_selenium_pourbaix(pH_steps=80, Eh_steps=80)

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


if __name__ == "__main__":
    test_write_selenium_pourbaix_outputs()
    print(f"Wrote Pourbaix plots to {POURBAIX_OUTPUT_DIR.resolve()}")
