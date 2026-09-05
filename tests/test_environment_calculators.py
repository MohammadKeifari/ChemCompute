import numpy as np
import pytest

from ChemCompute import Compound, Enviroment, EquilibriumResult, Reaction


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
    rxn = Reaction.from_string_simple_syntax(
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
