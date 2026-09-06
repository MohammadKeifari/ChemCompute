"""Tests for Enviroment: equilibrium, kinetics, composition, buffering, titration, activity, UV-Vis, bio."""

import math
import os

import numpy as np
import pytest

from ChemCompute import (
    ActivityModel,
    Compound,
    Enviroment,
    EquilibriumResult,
    HalfReaction,
    Pourbaix,
    Reaction,
    ScaledEnviroment,
    SpectrumSpec,
    Titration,
    XS,
    buffer_diagnostics,
    competitive_inhibition,
    ionic_strength,
    mix_sample_with_titrant,
    single_substrate_mm,
    uvvis_spectrum,
)

from helpers import ammonia_buffer_env, aq, simple_rxn, water_env, weak_acid_env


@pytest.fixture
def simple_equilibrium_environment():
    A = Compound("A")
    B = Compound("B")
    reactants = [{"stoichiometric_coefficient": 1, "compound": A, "rate_dependency": 1}]
    products = [{"stoichiometric_coefficient": 1, "compound": B, "rate_dependency": 1}]
    rxn = Reaction(reactants, products, [1.0], [0.0], K=2.0, kf=0.5, kb=0.25)
    return Enviroment(rxn, T=298)


@pytest.fixture
def simple_kinetic_environment():
    rxn = Reaction.from_string(
        "A > B",
        concentrations=[1.0, 0.0],
        K=2.0,
        kf=0.5,
        kb=0.25,
    )
    return Enviroment(rxn, T=298)


def test_env_equilibrium_default(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    result = env.equilibrium(max_iter=1000, tol=1e-8)
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)


@pytest.mark.parametrize("method", ["bgd", "sgd", "newton"])
def test_env_equilibrium_methods(simple_equilibrium_environment, method):
    env = simple_equilibrium_environment
    result = env.equilibrium(method=method, max_iter=1000, tol=1e-8)
    assert len(result) == 2
    assert np.isclose(result[0] + result[1], 1.0, atol=0.01)


@pytest.mark.parametrize("loss", ["log_quotient", "quotient_error", "log_huber"])
def test_env_equilibrium_loss_functions(simple_equilibrium_environment, loss):
    env = simple_equilibrium_environment
    result = env.equilibrium(method="bgd", loss=loss, max_iter=1000, tol=1e-6)
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)


def test_env_equilibrium_invalid_method(simple_equilibrium_environment):
    with pytest.raises(ValueError, match="Invalid method"):
        simple_equilibrium_environment.equilibrium(method="invalid")


def test_env_equilibrium_invalid_loss(simple_equilibrium_environment):
    with pytest.raises(ValueError, match="Invalid loss"):
        simple_equilibrium_environment.equilibrium(loss="invalid")


def test_compound_labels(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    assert env.compound_labels == ["A", "B"]
    assert env.compound_labels == [c.formula for c in env.compounds]


def test_concentrations_dict(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    assert env.concentrations_dict == {"A": 1.0, "B": 0.0}

    result = env.equilibrium(method="newton", tol=1e-10, return_details=True)
    assert result.concentrations_dict == dict(zip(result.compounds, result.concentrations))
    assert set(result.concentrations_dict) == {"A", "B"}


def test_equilibrium_return_details(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    result = env.equilibrium(method="newton", tol=1e-10, return_details=True)

    assert isinstance(result, EquilibriumResult)
    assert result.compounds == ["A", "B"]
    assert len(result.concentrations) == 2
    assert len(result.reaction_extents) == 1
    assert len(result.reaction_quotient_error) == 1
    assert len(result.reaction_quotient_ratio) == 1
    assert result.stop_reason in {"residual_tol", "quotient_error_limit", "max_iter"}
    assert result.iterations >= 1


def test_reaction_quotient_error_near_equilibrium(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    result = env.equilibrium(method="newton", tol=1e-12, return_details=True)

    assert result.max_reaction_quotient_error < 1e-4
    assert np.allclose(result.reaction_quotient_ratio, [1.0], rtol=1e-3, atol=1e-3)


def test_quotient_error_limit_stopping(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    result = env.equilibrium(
        method="bgd",
        quotient_error_limit=0.01,
        max_iter=5000,
        return_details=True,
    )

    assert result.stop_reason == "quotient_error_limit"
    assert result.max_reaction_quotient_error <= 0.01
    assert result.converged is True


def test_quotient_error_limit_ignores_tol(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    loose_tol_result = env.equilibrium(
        method="bgd",
        tol=1e30,
        quotient_error_limit=0.01,
        max_iter=5000,
        return_details=True,
    )
    assert loose_tol_result.stop_reason == "quotient_error_limit"
    assert loose_tol_result.max_reaction_quotient_error <= 0.01


def test_equilibrium_criterion_met_with_tol(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    result = env.equilibrium(method="newton", tol=1e-10, return_details=True)

    assert result.criterion_type == "residual_tol"
    assert result.criterion_met is True
    assert result.criterion_value < result.criterion_limit


def test_equilibrium_criterion_met_with_quotient_limit(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    result = env.equilibrium(
        method="bgd",
        quotient_error_limit=0.01,
        max_iter=5000,
        return_details=True,
    )

    assert result.criterion_type == "quotient_error"
    assert result.criterion_met is True
    assert result.criterion_value <= result.criterion_limit


def test_equilibrium_criterion_not_met_at_max_iter(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    result = env.equilibrium(method="bgd", tol=1e-30, max_iter=1, return_details=True)

    assert result.stop_reason == "max_iter"
    assert result.criterion_type == "residual_tol"
    assert result.criterion_met is False


def test_equilibrium_last_result_without_return_details(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    concentrations = env.equilibrium(method="newton", tol=1e-10)

    assert isinstance(concentrations, list)
    last = env.last_equilibrium_result
    assert isinstance(last, EquilibriumResult)
    assert last.concentrations == concentrations
    assert last.stop_reason in {"residual_tol", "quotient_error_limit", "max_iter"}
    assert last.iterations >= 1
    assert len(last.q_over_k) == 1
    assert np.isclose(last.q_over_k[0], 1.0, rtol=1e-3, atol=1e-3)


def test_apply_equilibrium_updates_concentrations(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    assert env.concentrations == [1.0, 0.0]

    result = env.apply_equilibrium(method="newton", tol=1e-10)

    assert isinstance(result, EquilibriumResult)
    assert env.concentrations == result.concentrations
    assert np.isclose(env.concentrations[0], 1.0 / 3.0, rtol=1e-4)
    assert np.isclose(env.concentrations[1], 2.0 / 3.0, rtol=1e-4)
    assert env.last_equilibrium_result is result
    assert result.iterations >= 1


def test_env_kinetics_basic(simple_kinetic_environment):
    env = simple_kinetic_environment
    result = env.kinetics(time=1.0, accuracy=1e-2, plot=False)
    assert len(result) >= 1
    final = result[-1]
    assert len(final) == 2
    assert all(c >= 0 for c in final)


def test_env_kinetics_mass_conservation(simple_kinetic_environment):
    env = simple_kinetic_environment
    initial_total = sum(env.concentrations)
    result = env.kinetics(time=2.0, accuracy=1e-2, plot=False)
    final_total = sum(result[-1])
    assert np.isclose(final_total, initial_total, atol=0.05)


def test_env_equilibrium_phase_exclusion():
    A_g = Compound("A", phase_point_list=[{"phase": "g", "temperature": 298}])
    B_l = Compound("B", phase_point_list=[{"phase": "l", "temperature": 298}])
    C_s = Compound("C", phase_point_list=[{"phase": "s", "temperature": 298}])
    rxn = Reaction(
        [{"stoichiometric_coefficient": 1, "compound": A_g, "rate_dependency": 1}],
        [
            {"stoichiometric_coefficient": 1, "compound": B_l, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": C_s, "rate_dependency": 1},
        ],
        [1.0],
        [0.0, 0.0],
        K=5.0,
        kf=0.3,
        kb=0.06,
    )
    env = Enviroment(rxn, T=298)
    result = env.equilibrium(max_iter=1000, tol=1e-8)
    assert len(result) == 3
    assert np.isclose(result[0], 0.2, atol=0.02)


@pytest.fixture
def weak_acid_trace_environment():
    ha = Compound("HA", phase_point_list=[{"phase": "aq", "temperature": 298}])
    hp = Compound("H+", phase_point_list=[{"phase": "aq", "temperature": 298}])
    am = Compound("A-", phase_point_list=[{"phase": "aq", "temperature": 298}])
    rxn = Reaction(
        [{"stoichiometric_coefficient": 1, "compound": ha, "rate_dependency": 1}],
        [
            {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": am, "rate_dependency": 1},
        ],
        [0.1],
        [1e-8, 1e-9],
        K=1e-5,
        kf=1.0,
        kb=1e5,
    )
    return Enviroment(rxn, T=298)


@pytest.fixture
def precipitation_solid_environment():
    ag = Compound("Ag+", phase_point_list=[{"phase": "aq", "temperature": 298}])
    cl = Compound("Cl-", phase_point_list=[{"phase": "aq", "temperature": 298}])
    agcl = Compound("AgCl", phase_point_list=[{"phase": "s", "temperature": 298}])
    rxn = Reaction(
        [
            {"stoichiometric_coefficient": 1, "compound": ag, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": cl, "rate_dependency": 1},
        ],
        [{"stoichiometric_coefficient": 1, "compound": agcl, "rate_dependency": 1}],
        [0.05, 0.05],
        [0.0],
        K=1e4,
        kf=1e4,
        kb=1.0,
    )
    return Enviroment(rxn, T=298)


@pytest.fixture
def solid_with_initial_amount_environment():
    ag = Compound("Ag+", phase_point_list=[{"phase": "aq", "temperature": 298}])
    cro4 = Compound("CrO4-2", phase_point_list=[{"phase": "aq", "temperature": 298}])
    precip = Compound("Ag2CrO4", phase_point_list=[{"phase": "s", "temperature": 298}])
    rxn = Reaction(
        [
            {"stoichiometric_coefficient": 2, "compound": ag, "rate_dependency": 2},
            {"stoichiometric_coefficient": 1, "compound": cro4, "rate_dependency": 1},
        ],
        [{"stoichiometric_coefficient": 1, "compound": precip, "rate_dependency": 1}],
        [0.2, 0.1],
        [0.5],
        K=1e6,
        kf=1e6,
        kb=1.0,
    )
    return Enviroment(rxn, T=298)


def test_weak_acid_equilibrium_with_trace_ions(weak_acid_trace_environment):
    env = weak_acid_trace_environment
    result = env.equilibrium(method="newton", tol=1e-12, return_details=True)

    assert result.concentrations[0] < 0.1
    assert result.concentrations[1] > 1e-8
    assert result.concentrations[2] > 1e-9
    assert np.isclose(result.concentrations[1], result.concentrations[2], rtol=1e-4)
    assert np.allclose(result.q_over_k, [1.0], rtol=1e-3, atol=1e-3)


def test_dilute_acid_with_1e_9_trace_concentrations():
    ha = Compound("HA", phase_point_list=[{"phase": "aq", "temperature": 298}])
    hp = Compound("H+", phase_point_list=[{"phase": "aq", "temperature": 298}])
    am = Compound("A-", phase_point_list=[{"phase": "aq", "temperature": 298}])
    rxn = Reaction(
        [{"stoichiometric_coefficient": 1, "compound": ha, "rate_dependency": 1}],
        [
            {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": am, "rate_dependency": 1},
        ],
        [1e-3],
        [1e-9, 1e-8],
        K=1e-4,
        kf=1.0,
        kb=1e4,
    )
    env = Enviroment(rxn, T=298)
    result = env.equilibrium(method="newton", tol=1e-12, return_details=True)

    assert all(c >= 0 for c in result.concentrations)
    assert result.concentrations[1] >= 1e-9
    assert result.max_reaction_quotient_error < 1e-6


def test_neutralization_liquid_water_excluded_from_quotient():
    hp = Compound("H+", phase_point_list=[{"phase": "aq", "temperature": 298}])
    ohm = Compound("OH-", phase_point_list=[{"phase": "aq", "temperature": 298}])
    h2o = Compound("H2O", phase_point_list=[{"phase": "l", "temperature": 298}])
    rxn = Reaction(
        [
            {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": ohm, "rate_dependency": 1},
        ],
        [{"stoichiometric_coefficient": 1, "compound": h2o, "rate_dependency": 1}],
        [1e-7, 1e-8],
        [55.0],
        K=1e14,
        kf=1e14,
        kb=1.0,
    )
    env = Enviroment(rxn, T=298)
    result = env.equilibrium(method="newton", tol=1e-12, return_details=True)

    assert np.allclose(result.q_over_k, [1.0], rtol=1e-3, atol=1e-3)
    assert result.concentrations[0] * result.concentrations[1] == pytest.approx(1e-14, rel=0.05)


def test_solid_and_liquid_excluded_from_quotient_matrix(precipitation_solid_environment):
    from ChemCompute._equilibrium import _build_context, _compute_lnQ

    env = precipitation_solid_environment
    ctx = _build_context(env, min_concentration=1e-12)

    assert ctx.A.shape == (1, 3)
    assert ctx.A[0, 0] == -1.0
    assert ctx.A[0, 1] == -1.0
    assert ctx.A[0, 2] == 0.0

    base_concentrations = np.array([0.01, 0.01, 0.04])
    lnQ_base = _compute_lnQ(ctx, base_concentrations)

    perturbed = base_concentrations.copy()
    perturbed[2] = 10.0
    lnQ_perturbed = _compute_lnQ(ctx, perturbed)

    assert np.isclose(lnQ_base, lnQ_perturbed)


def test_solid_product_accumulates_and_aqueous_depleted(precipitation_solid_environment):
    env = precipitation_solid_environment
    initial = env.concentrations.copy()
    result = env.equilibrium(method="newton", tol=1e-12, return_details=True)

    assert result.concentrations[0] < initial[0]
    assert result.concentrations[1] < initial[1]
    assert result.concentrations[2] > initial[2]
    assert np.allclose(result.q_over_k, [1.0], rtol=1e-3, atol=1e-3)


def test_pre_existing_solid_increases_during_precipitation(solid_with_initial_amount_environment):
    env = solid_with_initial_amount_environment
    initial_solid = env.concentrations[2]
    result = env.equilibrium(method="newton", tol=1e-12, return_details=True)

    assert result.concentrations[2] > initial_solid
    assert result.concentrations[0] < env.concentrations[0]
    assert np.allclose(result.q_over_k, [1.0], rtol=1e-3, atol=1e-3)


def test_excess_compound_concentration_fixed_during_equilibrium():
    from ChemCompute._equilibrium import _build_context

    h2o = Compound("H2O", phase_point_list=[{"phase": "l", "temperature": 298}], excess=True)
    hp = Compound("H+", phase_point_list=[{"phase": "aq", "temperature": 298}])
    ohm = Compound("OH-", phase_point_list=[{"phase": "aq", "temperature": 298}])
    rxn = Reaction(
        [{"stoichiometric_coefficient": 1, "compound": h2o, "rate_dependency": 1}],
        [
            {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": ohm, "rate_dependency": 1},
        ],
        [55.5],
        [1e-7, 1e-7],
        K=1e-14,
    )
    env = Enviroment(rxn, T=298)
    ctx = _build_context(env, min_concentration=1e-12)
    assert np.allclose(ctx.S[0, :], 0.0)

    result = env.equilibrium(method="newton", tol=1e-12, return_details=True)
    assert result.concentrations_dict["H2O"] == 55.5


def test_reaction_from_string_infinite_k():
    rxn = Reaction.from_string(
        "HCl.aq > H+ & Cl-",
        concentrations=[0.06, 0.0, 0.0],
        infinite_K=True,
    )
    assert rxn.infinite_K is True

    env = Enviroment(rxn, T=298)
    result = env.equilibrium(method="newton", tol=1e-12, return_details=True)
    assert result.concentrations_dict["HCl"] == pytest.approx(0.0, abs=1e-9)
    assert result.concentrations_dict["H+"] == pytest.approx(0.06, rel=1e-4)


def test_set_excess_with_xs_keeps_concentration():
    from ChemCompute import XS

    rxn = Reaction.from_string(
        "H2O.l > H+ & OH-",
        concentrations=[55.5, 1e-7, 1e-7],
        K=1e-14,
    )
    env = Enviroment(rxn, T=298)
    env.set_excess({"H2O": XS()})

    assert rxn.compounds[0]["excess"] is False
    assert env.compounds[0].excess is True
    assert env.concentrations_dict["H2O"] == 55.5

    result = env.equilibrium(method="newton", tol=1e-12, return_details=True)
    assert result.concentrations_dict["H2O"] == 55.5


def test_reaction_concentrations_dict_defaults_missing_to_zero():
    rxn = Reaction.from_string(
        "HA.aq > H+ & A-",
        concentrations={"HA": 0.1},
        K=1e-5,
    )
    assert rxn.compounds[0]["concentration"] == 0.1
    assert rxn.compounds[1]["concentration"] == 0.0
    assert rxn.compounds[2]["concentration"] == 0.0


def test_reaction_concentrations_dict_xs_marks_excess_on_entry():
    from ChemCompute import XS

    rxn = Reaction.from_string(
        "H2O.l > H+ & OH-",
        concentrations={"H2O": XS(55.5), "H+": 1e-7, "OH-": 1e-7},
        K=1e-14,
    )
    assert rxn.compounds[0]["excess"] is True
    assert rxn.compounds[0]["concentration"] == 55.5

    env = Enviroment(rxn, T=298)
    assert env.compounds[0].excess is True
    result = env.equilibrium(method="newton", tol=1e-12, return_details=True)
    assert result.concentrations_dict["H2O"] == 55.5


def test_reaction_concentrations_list_xs_marks_excess_on_entry():
    from ChemCompute import XS

    rxn = Reaction.from_string(
        "H2O.l > H+ & OH-",
        concentrations=[XS(55.5), 1e-7, 1e-7],
        K=1e-14,
    )
    assert rxn.compounds[0]["excess"] is True
    assert rxn.compounds[0]["concentration"] == 55.5

    env = Enviroment(rxn, T=298)
    assert env.compounds[0].excess is True
    result = env.equilibrium(method="newton", tol=1e-12, return_details=True)
    assert result.concentrations_dict["H2O"] == 55.5


def test_reaction_concentrations_list_xs_bare_marks_excess():
    from ChemCompute import XS

    rxn = Reaction.from_string(
        "CaF2.s > Ca+2 & 2_F-",
        concentrations=[XS(), 0.0, 0.0],
        K=5e-9,
    )
    assert rxn.compounds[0]["excess"] is True
    assert rxn.compounds[0]["concentration"] == 0.0


def test_concentrations_dict_accepts_xs_at_init():
    from ChemCompute import XS

    rxn = Reaction.from_string(
        "CaF2.s > Ca+2 & 2_F-",
        concentrations={"CaF2": XS(10.0)},
        K=5e-9,
    )
    env = Enviroment(rxn, T=298)
    assert rxn.compounds[0]["excess"] is True
    assert env.compounds[0].excess is True
    assert env.concentrations_dict["CaF2"] == 10.0


def test_infinite_k_drives_strong_acid_dissociation():
    hcl = Compound("HCl", phase_point_list=[{"phase": "aq", "temperature": 298}])
    hp = Compound("H+", phase_point_list=[{"phase": "aq", "temperature": 298}])
    clm = Compound("Cl-", phase_point_list=[{"phase": "aq", "temperature": 298}])
    rxn = Reaction(
        [{"stoichiometric_coefficient": 1, "compound": hcl, "rate_dependency": 1}],
        [
            {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": clm, "rate_dependency": 1},
        ],
        [0.06],
        [0.0, 0.0],
        infinite_K=True,
    )
    env = Enviroment(rxn, T=298)
    result = env.equilibrium(method="newton", tol=1e-12, return_details=True)

    assert result.concentrations_dict["HCl"] == pytest.approx(0.0, abs=1e-9)
    assert result.concentrations_dict["H+"] == pytest.approx(0.06, rel=1e-4)
    assert result.concentrations_dict["Cl-"] == pytest.approx(0.06, rel=1e-4)


def test_infinite_k_precipitation_consumes_limiting_ion():
    ag = Compound("Ag+", phase_point_list=[{"phase": "aq", "temperature": 298}])
    clm = Compound("Cl-", phase_point_list=[{"phase": "aq", "temperature": 298}])
    agcl = Compound("AgCl", phase_point_list=[{"phase": "s", "temperature": 298}])
    rxn = Reaction(
        [
            {"stoichiometric_coefficient": 1, "compound": ag, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": clm, "rate_dependency": 1},
        ],
        [{"stoichiometric_coefficient": 1, "compound": agcl, "rate_dependency": 1}],
        [0.01, 0.05],
        [0.0],
        infinite_K=True,
    )
    env = Enviroment(rxn, T=298)
    result = env.equilibrium(method="newton", tol=1e-12, return_details=True)

    assert result.concentrations_dict["Ag+"] == pytest.approx(0.0, abs=1e-12)
    assert result.concentrations_dict["AgCl"] == pytest.approx(0.01, rel=1e-4)
    assert result.concentrations_dict["Cl-"] == pytest.approx(0.04, rel=1e-4)


def test_fluoride_speciation_environment():
    import sys
    from pathlib import Path

    manual_dir = Path(__file__).resolve().parent / "manual"
    sys.path.insert(0, str(manual_dir))
    from complex_environments import COMPLEX_EQUILIBRIUM_KWARGS, ENV16_EXPECTED, build_env16

    env = build_env16()
    result = env.equilibrium(**COMPLEX_EQUILIBRIUM_KWARGS["env16"], return_details=True)

    assert result.criterion_met
    assert result.concentrations_dict["HF"] > 0.01
    assert result.concentrations_dict["Ca+2"] > 0.005
    ksp = result.concentrations_dict["Ca+2"] * result.concentrations_dict["F-"] ** 2
    assert ksp == pytest.approx(5e-9, rel=1e-4)
    assert result.concentrations_dict["HCl"] == pytest.approx(0.0, abs=1e-9)
    assert result.concentrations_dict["Cl-"] == pytest.approx(0.06, rel=1e-3)
    assert result.concentrations_dict["CaF2"] == 10.0
    assert result.concentrations_dict["H2O"] == 55.5
    for calc, exp in zip(result.concentrations, ENV16_EXPECTED):
        denom = max(abs(exp), abs(calc), 1e-12)
        assert abs(calc - exp) / denom <= 0.05


def test_coupled_network_sgd_from_zero_products():
    """SGD must converge when intermediate/product species start at zero."""
    a, b, c = Compound("A"), Compound("B"), Compound("C")
    env = Enviroment(
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": a, "rate_dependency": 1}],
            [{"stoichiometric_coefficient": 1, "compound": b, "rate_dependency": 1}],
            [1.0],
            [0.0],
            K=2.0,
            kf=0.5,
            kb=0.25,
        ),
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": b, "rate_dependency": 1}],
            [{"stoichiometric_coefficient": 1, "compound": c, "rate_dependency": 1}],
            [0.0],
            [0.0],
            K=1.5,
            kf=0.3,
            kb=0.2,
        ),
        T=298,
    )
    result = env.equilibrium(method="sgd", tol=1e-8, max_iter=8000, return_details=True)

    assert result.criterion_met
    assert np.allclose(result.concentrations, [1.0 / 6.0, 1.0 / 3.0, 0.5], rtol=0.03)


def test_water_autoionization_uses_kw_not_bulk_molarity():
    """Excess liquid water uses Kw = [H+][OH-], not K * [H2O]."""
    import sys
    from pathlib import Path

    manual_dir = Path(__file__).resolve().parent / "manual"
    sys.path.insert(0, str(manual_dir))
    from complex_environments import COMPLEX_EQUILIBRIUM_KWARGS, build_env16

    result = build_env16().equilibrium(
        **COMPLEX_EQUILIBRIUM_KWARGS["env16"], return_details=True
    )
    kw = result.concentrations_dict["H+"] * result.concentrations_dict["OH-"]
    assert kw == pytest.approx(1e-14, rel=1e-6)
    water_index = 3
    assert result.reaction_quotient_ratio[water_index] == pytest.approx(1.0, rel=1e-6)
    assert result.reaction_quotient_error[water_index] == pytest.approx(0.0, abs=1e-9)


# --- Composition ---


class TestFromCompounds:
    def test_compounds_only_no_reactions(self):
        env = Enviroment.from_compounds({"Na+": 0.1, "Cl-": 0.1}, volume=1.0)
        assert env.reactions == []
        assert np.isclose(env.concentrations_dict["Na+"], 0.1)
        assert np.isclose(env.concentrations_dict["Cl-"], 0.1)
        assert env.volume == 1.0

    def test_equilibrium_no_reactions(self):
        env = Enviroment.from_compounds({"A": 0.5}, volume=1.0)
        result = env.equilibrium(return_details=True)
        assert result.stop_reason == "no_reactions"
        assert result.concentrations == [0.5]


class TestConcentrationOverride:
    def test_override_beats_reaction_sum(self):
        rxn = simple_rxn(initial_a=1.0, initial_b=0.0)
        env = Enviroment(rxn, concentrations={"A": 0.25})
        assert np.isclose(env.concentrations_dict["A"], 0.25)

    def test_override_adds_new_species(self):
        rxn = simple_rxn()
        extra = aq("X")
        env = Enviroment(rxn, concentrations={extra: 0.03})
        assert "X" in env.concentrations_dict
        assert np.isclose(env.concentrations_dict["X"], 0.03)


class TestStandaloneReaction:
    def test_equilibrium_matches_env_wrapper(self):
        rxn = simple_rxn()
        direct = rxn.equilibrium(method="newton", tol=1e-10)
        wrapped = Enviroment(rxn).equilibrium(method="newton", tol=1e-10)
        assert np.allclose(direct, wrapped)

    def test_kinetics_runs(self):
        rxn = Reaction.from_string(
            "A > B",
            [1.0, 0.0],
            K=2.0,
            kf=0.5,
            kb=0.25,
        )
        checkpoints = rxn.kinetics(time=1.0, accuracy=0.1, plot=False)
        assert checkpoints[-1][0] < 1.0


class TestEnvironmentMixing:
    def test_add_equal_volumes(self):
        envA = Enviroment.from_compounds({"A": 1.0}, volume=1.0)
        envB = Enviroment.from_compounds({"B": 0.1}, volume=1.0)
        envC = envA + envB
        assert envC.volume == 2.0
        assert np.isclose(envC.concentrations_dict["A"], 0.5)
        assert np.isclose(envC.concentrations_dict["B"], 0.05)

    def test_coefficient_mix(self):
        envA = Enviroment.from_compounds({"A": 1.0}, volume=1.0)
        envB = Enviroment.from_compounds({"B": 0.1}, volume=1.0)
        envD = 0.5 * envA + 4 * envB
        assert np.isclose(envD.volume, 4.5)
        assert np.isclose(envD.concentrations_dict["A"], 0.5 / 4.5)
        assert np.isclose(envD.concentrations_dict["B"], 0.4 / 4.5)

    def test_combine_matches_operator(self):
        envA = Enviroment.from_compounds({"A": 1.0}, volume=1.0)
        envB = Enviroment.from_compounds({"B": 0.1}, volume=1.0)
        env_op = 0.5 * envA + 4 * envB
        env_fn = Enviroment.combine((0.5, envA), (4.0, envB))
        assert np.isclose(env_op.volume, env_fn.volume)
        assert env_op.concentrations_dict == env_fn.concentrations_dict

    def test_add_compounds_slug(self):
        envB = Enviroment.from_compounds({"B": 0.1}, volume=1.0)
        envC = envB.add_compounds({"A": 1.0}, volume=1.0)
        assert envC.volume == 2.0
        assert np.isclose(envC.concentrations_dict["A"], 0.5)
        assert np.isclose(envC.concentrations_dict["B"], 0.05)

    def test_add_compounds_with_coefficient(self):
        envB = Enviroment.from_compounds({"B": 0.1}, volume=1.0)
        envE = envB.add_compounds({"A": 1.0}, volume=1.0, coefficient=4.0)
        assert np.isclose(envE.volume, 5.0)

    def test_temperature_mismatch_raises(self):
        envA = Enviroment.from_compounds({"A": 1.0}, T=298)
        envB = Enviroment.from_compounds({"B": 0.1}, T=310)
        with pytest.raises(ValueError):
            _ = envA + envB

    def test_invalid_coefficient_raises(self):
        envA = Enviroment.from_compounds({"A": 1.0})
        with pytest.raises(ValueError):
            _ = 0 * envA + envA

    def test_scaled_type(self):
        envA = Enviroment.from_compounds({"A": 1.0})
        assert isinstance(0.5 * envA, ScaledEnviroment)


# --- Buffering (enforcement) ---


class TestConstantPHEquilibrium:
    def test_buffered_h_plus_unchanged(self):
        env = weak_acid_env(h_plus=1e-7, buffer=["H+"])
        result = env.equilibrium(return_details=True)
        assert np.isclose(result.concentrations_dict["H+"], 1e-7, rtol=1e-3)
        assert result.concentrations_dict["HA"] < 0.1
        assert result.concentrations_dict["A-"] > 0.0

    def test_unbuffered_h_plus_shifts(self):
        env = weak_acid_env(h_plus=1e-7, buffer=None)
        result = env.equilibrium(return_details=True)
        assert result.concentrations_dict["H+"] > 1e-7

    def test_q_over_k_finite_with_buffered_h_plus(self):
        env = weak_acid_env(h_plus=1e-7, buffer=["H+"])
        result = env.equilibrium(return_details=True)
        assert all(math.isfinite(q) for q in result.q_over_k)

    def test_explicit_buffer_target(self):
        env = weak_acid_env(h_plus=1e-7, buffer={"H+": 1e-5})
        assert np.isclose(env.concentrations_dict["H+"], 1e-5)
        result = env.equilibrium(return_details=True)
        assert np.isclose(result.concentrations_dict["H+"], 1e-5, rtol=1e-3)

    def test_unknown_buffered_species_raises(self):
        env = weak_acid_env(buffer=None)
        with pytest.raises(ValueError, match="Buffered species"):
            env.set_buffer(["MissingIon"])


class TestConstantPHKinetics:
    def test_h_plus_flat_during_integration(self):
        env = weak_acid_env(h_plus=1e-7, buffer=["H+"])
        traces = env.kinetics(time=1.0, accuracy=0.1)
        h_index = env.compound_labels.index("H+")
        for snapshot in traces:
            assert np.isclose(snapshot[h_index], 1e-7, rtol=1e-6)


class TestBufferCombine:
    def test_merged_h_plus_target_from_mix(self):
        env_a = Enviroment.from_compounds(
            {"H+": 1e-7, "Cl-": 0.1},
            volume=1.0,
            buffer=["H+"],
        )
        env_b = Enviroment.from_compounds(
            {"H+": 1e-5, "Na+": 0.1},
            volume=1.0,
            buffer=["H+"],
        )
        combined = env_a + env_b
        expected = (1e-7 + 1e-5) / 2.0
        assert np.isclose(combined.buffer_targets["H+"], expected)
        assert np.isclose(combined.concentrations_dict["H+"], expected)


class TestSetBuffer:
    def test_set_buffer_resnapshots_h_plus(self):
        env = Enviroment.from_compounds({"H+": 1e-7, "Cl-": 0.1})
        updated = env.concentrations
        updated[env.compound_labels.index("H+")] = 1e-5
        env.concentrations = updated
        env.set_buffer(["H+"])
        assert np.isclose(env.buffer_targets["H+"], 1e-5)

    def test_copy_preserves_buffer_targets(self):
        env = Enviroment.from_compounds({"H+": 1e-7, "Cl-": 0.1}, buffer=["H+"])
        copied = env.copy()
        assert copied.buffer_targets == env.buffer_targets


class TestGenericBufferMechanism:
    def test_buffered_species_fixed_in_simple_equilibrium(self):
        a = aq("A")
        b = aq("B")
        rxn = Reaction(
            reactants=[{"stoichiometric_coefficient": 1, "compound": a, "rate_dependency": 1}],
            products=[{"stoichiometric_coefficient": 1, "compound": b, "rate_dependency": 1}],
            reactants_concentration=[0.5],
            products_concentration=[0.0],
            K=2.0,
        )
        env = Enviroment(rxn, buffer=["A"])
        result = env.equilibrium(return_details=True)
        assert np.isclose(result.concentrations_dict["A"], 0.5)
        assert result.concentrations_dict["B"] > 0.0


# --- Titration ---


class TestTitration:
    def test_strong_base_titration_pH_rises(self):
        sample = water_env(h_plus=1e-3, oh_minus=1e-11, volume=0.1)
        titrant = Enviroment.from_compounds({"OH-": 0.1, "Na+": 0.1}, volume=1.0)
        result = Titration(
            sample,
            titrant,
            volume_min=0.0,
            volume_max=0.02,
            steps=20,
        ).run(method="newton", tol=1e-8)

        assert len(result.titrant_volumes) == 20
        assert result.pH[-1] > result.pH[0]
        assert result.matrix().shape == (20, len(result.compound_labels))

    def test_matrix_rows_match_volumes(self):
        sample = water_env(volume=0.1)
        titrant = Enviroment.from_compounds({"OH-": 0.01, "Na+": 0.01}, volume=1.0)
        volumes = [0.0, 0.005, 0.01]
        result = Titration(sample, titrant, volumes=volumes).run(method="newton", tol=1e-8)

        assert result.titrant_volumes == volumes
        assert result.total_volumes[0] == sample.volume
        assert result.total_volumes[-1] > sample.volume

    def test_species_series(self):
        sample = water_env(volume=0.1)
        titrant = Enviroment.from_compounds({"OH-": 0.01, "Na+": 0.01}, volume=1.0)
        result = Titration(sample, titrant, volume_min=0.0, volume_max=0.01, steps=5).run(
            method="newton",
            tol=1e-8,
        )

        h_series = result.species("H+")
        assert len(h_series) == 5
        assert np.all(np.diff(h_series) <= 0)

    def test_species_unknown_raises(self):
        sample = water_env(volume=0.1)
        titrant = Enviroment.from_compounds({"OH-": 0.01}, volume=1.0)
        result = Titration(sample, titrant, volume_min=0.0, volume_max=0.001, steps=3).run(
            method="newton",
            tol=1e-8,
        )
        with pytest.raises(KeyError):
            result.species("NotASpecies")

    def test_sample_not_mutated(self):
        sample = water_env(h_plus=1e-3, volume=0.1)
        initial = sample.concentrations_dict.copy()
        titrant = Enviroment.from_compounds({"OH-": 0.1, "Na+": 0.1}, volume=1.0)
        Titration(sample, titrant, volume_min=0.0, volume_max=0.01, steps=5).run(
            method="newton",
            tol=1e-8,
        )
        assert sample.concentrations_dict == initial
        assert sample.volume == 0.1

    def test_mix_zero_titrant_volume_returns_sample_copy(self):
        sample = water_env(volume=0.1)
        titrant = Enviroment.from_compounds({"OH-": 0.1}, volume=1.0)
        mixed = mix_sample_with_titrant(sample, titrant, 0.0)
        assert mixed.volume == sample.volume
        assert mixed.concentrations_dict == sample.concentrations_dict
        assert mixed is not sample

    def test_temperature_mismatch_raises(self):
        sample = water_env()
        titrant = Enviroment.from_compounds({"OH-": 0.1}, T=310, volume=1.0)
        with pytest.raises(ValueError, match="temperature"):
            Titration(sample, titrant)

    def test_plot_runs_without_display(self):
        sample = water_env(volume=0.1)
        titrant = Enviroment.from_compounds({"OH-": 0.01, "Na+": 0.01}, volume=1.0)
        result = Titration(sample, titrant, volume_min=0.0, volume_max=0.01, steps=5).run(
            method="newton",
            tol=1e-8,
        )
        result.plot(species=["H+", "OH-"], plot="save", directory="titration_test_plot.png")
        result.plot_pH(plot="save", directory="titration_test_pH.png")
        assert os.path.isfile("titration_test_plot.png")
        assert os.path.isfile("titration_test_pH.png")
        os.remove("titration_test_plot.png")
        os.remove("titration_test_pH.png")


class TestAmmoniumTartrateSilverFluorideTitration:
    """50 mL of 0.07 M AgF titrated with 0.01 M (NH4)2T until Ag2T would precipitate."""

    def test_qsp_surpasses_ksp_near_0_004_microliters(self):
        ksp = 4e-12
        eq_kwargs = dict(method="newton", tol=1e-10, max_iter=8000, min_concentration=1e-20)

        agf = Enviroment(
            Reaction.from_string(
                "H2O.l > H+ & OH-",
                concentrations={"H2O": XS(55.5), "H+": 1e-7, "OH-": 1e-7},
                K=1e-14,
            ),
            Reaction.from_string("HF.aq > H+ & F-", K=10 ** (-3.1)),
            concentrations={"Ag+": 0.07, "F-": 0.07},
            volume=0.050,
        )
        titrant = Enviroment(
            Reaction.from_string("H2T > H+ & HT-", K=10 ** (-4.2)),
            Reaction.from_string("HT- > H+ & T-2", K=10 ** (-6.5)),
            Reaction.from_string("NH4+ > H+ & NH3", K=10 ** (-9.2)),
            Reaction.from_string("Ag+ & 2_NH3 > Ag(NH3)2+", K=2e7),
            concentrations={"NH4+": 0.02, "T-2": 0.01},
            volume=1.0,
        )

        last_safe_ul = None
        first_precip_ul = None
        previous_qsp = None
        for volume_ul in (0.001, 0.003, 0.004, 0.0042, 0.0045, 0.01):
            mixed = agf + (volume_ul * 1e-6) * titrant
            result = mixed.equilibrium(**eq_kwargs, return_details=True)
            assert result.criterion_met
            qsp = result.concentrations_dict["Ag+"] ** 2 * result.concentrations_dict["T-2"]
            if previous_qsp is not None:
                assert qsp > previous_qsp
            previous_qsp = qsp
            if qsp <= ksp:
                last_safe_ul = volume_ul
            elif first_precip_ul is None:
                first_precip_ul = volume_ul

        assert last_safe_ul == 0.0042
        assert first_precip_ul == 0.0045


class TestFeSCNUnknownSolutions:
    """Identify A–D from 470 nm absorbances, then predict the equal-volume mix."""

    def test_table_identifies_abcd_and_equal_volume_absorbance(self):
        wavelength = 470e-9
        eq_kwargs = dict(method="newton", tol=1e-10, max_iter=8000, min_concentration=1e-20)

        fe = Enviroment(
            Reaction.from_string("Fe+3 > FeOH+2 & H+", K=10 ** (-2.90)),
            Reaction.from_string("Fe+3 & SCN- > FeSCN+2", K=10 ** 3.22),
            concentrations={"Fe+3": 0.001, "NO3-": 0.003},
            volume=1.0,
        )
        fe.set_spectrum("FeSCN+2", SpectrumSpec(points=[(wavelength, 3700.0)], extrapolate="none"))
        scn = Enviroment.from_compounds({"K+": 0.0005, "SCN-": 0.0005}, volume=1.0)
        acid = Enviroment(
            Reaction.from_string(
                "HNO3.aq > H+ & NO3-",
                concentrations={"HNO3": 0.120},
                infinite_K=True,
            ),
            volume=1.0,
        )
        water = Enviroment(
            Reaction.from_string(
                "H2O.l > H+ & OH-",
                concentrations={"H2O": XS(55.5), "H+": 1e-7, "OH-": 1e-7},
                K=1e-14,
            ),
            volume=1.0,
        )

        solution_a = fe
        solution_b = water
        solution_c = acid
        solution_d = scn

        cases = (
            (20, 20, 40, 20, 0.082, 0.005),
            (20, 55, 5, 20, 0.073, 0.005),
            (10, 40, 40, 10, 0.024, 0.005),
            (50, 10, 5, 35, 0.242, 0.005),
            (20, 10, 10, 50, 0.173, 0.04),
        )
        for va, vb, vc, vd, measured, abs_tol in cases:
            mixed = (
                (va * 1e-3) * solution_a
                + (vb * 1e-3) * solution_b
                + (vc * 1e-3) * solution_c
                + (vd * 1e-3) * solution_d
            )
            result = mixed.apply_equilibrium(**eq_kwargs)
            assert result.criterion_met
            absorbance = uvvis_spectrum(
                mixed,
                wavelengths=[wavelength],
                path_length=1.0,
            )[0]
            assert absorbance == pytest.approx(measured, abs=abs_tol)

        equal = (
            0.025 * solution_a
            + 0.025 * solution_b
            + 0.025 * solution_c
            + 0.025 * solution_d
        )
        result = equal.apply_equilibrium(**eq_kwargs)
        assert result.criterion_met
        absorbance = uvvis_spectrum(equal, wavelengths=[wavelength], path_length=1.0)[0]
        assert absorbance == pytest.approx(0.119, abs=0.002)


class TestStrontiumFluoridePurification:
    """5.00 g of 90% SrF2 / 10% PbF2 leached in acetic acid or in KI."""

    def test_remaining_solid_mass_and_srf2_purity_for_both_methods(self):
        mw_srf2 = 87.6 + 2 * 19.0
        mw_pbf2 = 207.2 + 2 * 19.0
        mw_pbi2 = 207.2 + 2 * 126.9
        n_srf2 = 0.90 * 5.00 / mw_srf2
        n_pbf2 = 0.10 * 5.00 / mw_pbf2
        v_leach = 0.200
        v_sample = 1e-9
        eq_kwargs = dict(method="newton", tol=1e-10, max_iter=200, min_concentration=1e-20)

        sample = Enviroment(
            Reaction.from_string(
                "SrF2.s > Sr+2 & 2_F-",
                concentrations={"SrF2": n_srf2 / v_sample},
                K=4e-9,
            ),
            Reaction.from_string(
                "PbF2.s > Pb+2 & 2_F-",
                concentrations={"PbF2": n_pbf2 / v_sample},
                K=3e-8,
            ),
            volume=v_sample,
        )
        acetic = Enviroment(
            Reaction.from_string(
                "H2O.l > H+ & OH-",
                concentrations={"H2O": XS(55.5), "H+": 1e-7, "OH-": 1e-7},
                K=1e-14,
            ),
            Reaction.from_string(
                "HOAc > H+ & OAc-",
                concentrations={"HOAc": 0.1},
                K=10 ** (-4.75),
            ),
            Reaction.from_string("HF.aq > H+ & F-", K=10 ** (-3.10)),
            concentrations={"Sr+2": 1e-6, "Pb+2": 1e-6, "F-": 1e-6, "HF": 1e-6, "OAc-": 1e-6},
            volume=v_leach,
        )
        ki = Enviroment(
            Reaction.from_string(
                "H2O.l > H+ & OH-",
                concentrations={"H2O": XS(55.5), "H+": 1e-7, "OH-": 1e-7},
                K=1e-14,
            ),
            Reaction.from_string("HF.aq > H+ & F-", K=10 ** (-3.10)),
            Reaction.from_string("PbI2.s > Pb+2 & 2_I-", K=9e-9),
            Reaction.from_string("Pb+2 & 4_I- > PbI4-2", K=10 ** 4.6),
            concentrations={
                "K+": 0.06,
                "I-": 0.06,
                "Sr+2": 1e-6,
                "Pb+2": 1e-6,
                "F-": 1e-6,
                "HF": 1e-6,
                "PbI4-2": 1e-6,
            },
            volume=v_leach,
        )

        def remaining_solids(env):
            conc = env.concentrations_dict
            volume = env.volume
            masses = {
                "SrF2": max(conc.get("SrF2", 0.0), 0.0) * volume * mw_srf2,
                "PbF2": max(conc.get("PbF2", 0.0), 0.0) * volume * mw_pbf2,
                "PbI2": max(conc.get("PbI2", 0.0), 0.0) * volume * mw_pbi2,
            }
            total = sum(masses.values())
            purity = 100.0 * masses["SrF2"] / total
            return masses, total, purity

        method1 = acetic + sample
        result1 = method1.apply_equilibrium(**eq_kwargs)
        assert result1.criterion_met
        masses1, total1, purity1 = remaining_solids(method1)
        assert total1 == pytest.approx(4.865, abs=0.005)
        assert purity1 == pytest.approx(92.33, abs=0.05)
        assert masses1["PbI2"] == pytest.approx(0.0, abs=1e-9)
        assert masses1["SrF2"] > masses1["PbF2"]

        method2 = ki + sample
        method2.apply_equilibrium(**eq_kwargs)
        masses2, total2, purity2 = remaining_solids(method2)
        assert masses2["PbF2"] == pytest.approx(0.0, abs=1e-6)
        assert masses2["PbI2"] == pytest.approx(0.940, abs=0.005)
        assert masses2["SrF2"] == pytest.approx(4.500, abs=0.005)
        assert total2 == pytest.approx(5.440, abs=0.005)
        assert purity2 == pytest.approx(82.73, abs=0.05)


class TestRedoxMixtureVsSCE:
    """Mix A3+/A+, D4+/D2+, and B3+/B+ solutions; closed-system potential vs SCE."""

    def test_mixture_potential_vs_sce_and_final_b3(self):
        from ChemCompute._half_reaction import F_FARADAY, R_GAS

        n_electrons = 2.0
        k_ad = math.exp(n_electrons * F_FARADAY * (1.5206 - 1.511) / (R_GAS * 298.0))
        k_ab = math.exp(n_electrons * F_FARADAY * (1.5206 - 1.5931) / (R_GAS * 298.0))
        eq_kwargs = dict(method="newton", tol=1e-12, max_iter=8000, min_concentration=1e-20)

        redox = Enviroment(
            Reaction.from_string("A+3 & D+2 > A+ & D+4", K=k_ad),
            Reaction.from_string("A+3 & B+ > A+ & B+3", K=k_ab),
            volume=1e-9,
        )
        mixed = (
            Enviroment.from_compounds({"A+3": 5e-3}, volume=0.050)
            + Enviroment.from_compounds({"D+2": 2e-3}, volume=0.050)
            + Enviroment.from_compounds({"A+": 3.5e-3}, volume=0.050)
            + Enviroment.from_compounds({"B+": 1.5e-3}, volume=0.100)
            + redox
        )
        result = mixed.apply_equilibrium(**eq_kwargs)
        assert result.criterion_met
        conc = result.concentrations_dict
        eh_she = 1.5206 + (R_GAS * mixed.T) / (n_electrons * F_FARADAY) * math.log(
            conc["A+3"] / conc["A+"]
        )
        eh_d = 1.511 + (R_GAS * mixed.T) / (n_electrons * F_FARADAY) * math.log(
            conc["D+4"] / conc["D+2"]
        )
        eh_b = 1.5931 + (R_GAS * mixed.T) / (n_electrons * F_FARADAY) * math.log(
            conc["B+3"] / conc["B+"]
        )
        assert eh_she == pytest.approx(eh_d, abs=1e-8)
        assert eh_she == pytest.approx(eh_b, abs=1e-8)
        assert eh_she - 0.241 == pytest.approx(1.2765, abs=1e-4)
        assert conc["B+3"] == pytest.approx(1.662e-6, rel=1e-3)
        assert mixed.volume == pytest.approx(0.250, abs=1e-8)


# --- Activity ---


class TestActivityModel:
    def test_davies_gamma_at_ionic_strength(self):
        env = Enviroment.__new__(Enviroment)
        na = aq("Na+", charge=1)
        cl = aq("Cl-", charge=-1)
        env.compounds = [na, cl]
        env.charge_map = {}
        env._T = 298
        conc = np.array([0.1, 0.1])
        model = ActivityModel("davies")
        gammas = model.gamma_array(env, conc, 298)
        assert gammas[0] < 1.0
        assert gammas[1] < 1.0
        assert np.isclose(gammas[0], gammas[1], rtol=0.05)

    def test_ionic_strength_na_cl(self):
        env = Enviroment.__new__(Enviroment)
        env.compounds = [aq("Na+", charge=1), aq("Cl-", charge=-1)]
        env.charge_map = {}
        env._T = 298
        i = ionic_strength(env, np.array([0.2, 0.2]), 298)
        assert np.isclose(i, 0.2)

    def test_equilibrium_with_davies_activity_model(self, simple_equilibrium_environment):
        env = simple_equilibrium_environment
        env.activity_model = "davies"
        env.charge_map = {"A": 0, "B": 0}
        result = env.equilibrium(method="newton", tol=1e-8, return_details=True)
        assert result.criterion_met
        assert all(c >= 0 for c in result.concentrations)


# --- Buffer diagnostics ---


class TestBufferDiagnostics:
    def test_buffer_beta_positive(self):
        env = ammonia_buffer_env()
        result = env.equilibrium(method="newton", tol=1e-10, return_details=True)
        diag = buffer_diagnostics(env, result.concentrations)
        assert diag.beta > 0
        assert not math.isnan(diag.pH)


# --- UV-Vis ---


class TestUVVis:
    def test_piecewise_flat_extrapolation(self):
        spec = SpectrumSpec(
            points=[(400e-9, 1000.0), (500e-9, 2000.0)],
            extrapolate="flat",
        )
        assert spec.epsilon(350e-9) == 1000.0
        assert spec.epsilon(450e-9) == 1500.0
        assert spec.epsilon(600e-9) == 2000.0

    def test_extrapolate_none(self):
        spec = SpectrumSpec(points=[(450e-9, 25000.0)], extrapolate="none")
        assert spec.epsilon(450e-9) == 25000.0
        assert spec.epsilon(500e-9) == 0.0

    def test_beer_lambert_linear(self):
        dye = aq("D")
        dummy = aq("X")
        rxn = Reaction(
            reactants=[{"stoichiometric_coefficient": 1, "compound": dye, "rate_dependency": 1}],
            products=[{"stoichiometric_coefficient": 1, "compound": dummy, "rate_dependency": 1}],
            reactants_concentration=[0.001],
            products_concentration=[0.0],
            K=1.0,
        )
        env = Enviroment(rxn)
        env.concentrations = [0.001, 0.0]
        env.set_spectrum("D", SpectrumSpec(points=[(500e-9, 1000.0)], extrapolate="flat"))
        a1 = uvvis_spectrum(env, wavelengths=[500e-9], path_length=0.01)[0]
        env.concentrations = [0.002, 0.0]
        a2 = uvvis_spectrum(env, wavelengths=[500e-9], path_length=0.01)[0]
        assert np.isclose(a2, 2 * a1)


# --- Phase, half-reactions, Pourbaix ---


def test_undetermined_phase_stays_in_quotient():
    from ChemCompute._equilibrium import _build_context, _compute_lnQ

    a = Compound("A")
    b = Compound("B")
    rxn = Reaction(
        [{"stoichiometric_coefficient": 1, "compound": a, "rate_dependency": 1}],
        [{"stoichiometric_coefficient": 1, "compound": b, "rate_dependency": 1}],
        [0.5],
        [0.5],
        K=1.0,
    )
    env = Enviroment(rxn)
    ctx = _build_context(env, min_concentration=1e-12)
    assert ctx.A[0, 0] != 0.0
    assert ctx.A[0, 1] != 0.0
    base = np.array([0.5, 0.5])
    perturbed = np.array([0.5, 2.0])
    assert not np.isclose(_compute_lnQ(ctx, base), _compute_lnQ(ctx, perturbed))


def test_half_reaction_string_concentrations_and_e_at():
    from ChemCompute import HalfReaction

    hr = HalfReaction.from_string(
        "Fe+3 & @e = Fe+2",
        concentrations=[0.01, 0.001],
        E0=0.771,
    )
    assert hr.n_electrons == 1.0
    assert hr.oxidized_concentration == [0.01]
    assert hr.reduced_concentration == [0.001]
    assert hr.E_at() == pytest.approx(0.83013, rel=1e-3)


def test_half_reaction_complex_syntax_and_h_plus_slope():
    from ChemCompute import HalfReaction

    hr = HalfReaction.from_string(
        "Fe(OH)3.s & 3_H+ & @e = Fe+2 & 3_H2O.l",
        concentrations=[1.0, 1e-7, 0.05, 1.0],
        E0=-0.55,
    )
    assert hr.net_h_plus_stoichiometry() == 3.0
    h_plus = next(entry for entry in hr.oxidized if entry["compound"].formula == "H+")
    assert h_plus["stoichiometric_coefficient"] == 3
    assert h_plus["rate_dependency"] == 3
    e_low = hr.E_at_pH(0.0)
    e_high = hr.E_at_pH(7.0)
    assert e_high < e_low


def test_reaction_rejects_electron_formula():
    with pytest.raises(ValueError, match="Electrons belong"):
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": Compound("e-"), "rate_dependency": 1}],
            [{"stoichiometric_coefficient": 1, "compound": Compound("A"), "rate_dependency": 1}],
            [1.0],
            [0.0],
        )


def test_env_mixed_reaction_and_half_reaction_constructor():
    from ChemCompute import HalfReaction

    hr = HalfReaction.from_string(
        "Fe+3 & @e = Fe+2",
        concentrations=[0.01, 0.001],
        E0=0.771,
    )
    env = Enviroment(hr)
    assert len(env.half_reactions) == 1
    assert env.compound_labels == ["Fe+3", "Fe+2"]


def test_fixed_electrode_potential_equilibrium():
    from ChemCompute import HalfReaction

    hr = HalfReaction.from_string(
        "Fe+3 & @e = Fe+2",
        concentrations=[0.01, 0.001],
        E0=0.771,
    )
    env = Enviroment(hr)
    env.set_electrode_potential(Eh=0.44)
    result = env.equilibrium(method="newton", return_details=True)
    assert result.electrode_Eh == pytest.approx(0.44)
    assert np.allclose(result.q_over_k, [1.0], rtol=1e-3, atol=1e-3)


def test_coupled_electrode_potential_two_half_reactions():
    from ChemCompute import HalfReaction

    hr_fe = HalfReaction.from_string(
        "Fe+3 & @e = Fe+2",
        concentrations=[0.01, 0.001],
        E0=0.771,
    )
    hr_ce = HalfReaction.from_string(
        "Ce+4 & @e = Ce+3",
        concentrations=[0.01, 0.001],
        E0=1.72,
    )
    env = Enviroment(hr_fe, hr_ce)
    result = env.equilibrium(method="newton", return_details=True)
    assert result.electrode_Eh is not None
    assert np.allclose(result.q_over_k, [1.0, 1.0], rtol=1e-3, atol=1e-3)


def test_duplicate_half_reaction_registration_raises():
    from ChemCompute import HalfReaction

    hr1 = HalfReaction.from_string("Fe+3 & @e = Fe+2", E0=0.771)
    hr2 = HalfReaction.from_string("Fe+3 & @e = Fe+2", E0=0.5)
    env = Enviroment()
    env.register_half_reactions([hr1])
    with pytest.raises(ValueError, match="Duplicate"):
        env.register_half_reactions([hr2])


def test_kinetics_warns_when_half_reactions_present():
    from ChemCompute import HalfReaction
    import warnings

    hr = HalfReaction.from_string("Fe+3 & @e = Fe+2", E0=0.771)
    env = Enviroment(
        Reaction.from_string("A > B", concentrations=[1.0, 0.0], K=1.0),
        hr,
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        env.kinetics(time=0.01, checkpoint_time=[0.01], plot=False)
    assert any("Half-reactions" in str(w.message) for w in caught)


def test_copy_preserves_half_reactions_and_electrode_eh():
    from ChemCompute import HalfReaction

    hr = HalfReaction.from_string(
        "Fe+3 & @e = Fe+2",
        concentrations=[0.01, 0.001],
        E0=0.771,
    )
    env = Enviroment(hr, electrode_Eh=0.25)
    copied = env.copy()
    assert copied.electrode_Eh == 0.25
    assert len(copied.half_reactions) == 1
    assert copied.half_reactions[0]._reaction_index is not None


def test_pourbaix_grid_run():
    from ChemCompute import HalfReaction, Pourbaix

    hr = HalfReaction.from_string(
        "Fe+3 & @e = Fe+2",
        concentrations=[0.01, 0.001],
        E0=0.771,
    )
    water = Compound("H2O", excess=True)
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
    assert len(diagram.boundary_lines) == 1
    assert diagram.speciation is not None
    assert diagram.track_species == ["Fe+3", "Fe+2"]


# --- Bio kinetics ---


class TestBioKinetics:
    def test_mm_template_runs(self):
        env = single_substrate_mm(s0=1.0, Vmax=1e-5, Km=1e-4)
        checkpoints = env.kinetics(time=10.0, accuracy=0.1, plot=False)
        final = checkpoints[-1]
        assert final[0] < 1.0
        assert final[1] > 0.0

    def test_competitive_template(self):
        env = competitive_inhibition(s0=1.0, i0=0.5, Ki=1e-4, Vmax=1e-5, Km=1e-4)
        checkpoints = env.kinetics(time=5.0, accuracy=0.1, plot=False)
        assert checkpoints[-1][0] < 1.0
