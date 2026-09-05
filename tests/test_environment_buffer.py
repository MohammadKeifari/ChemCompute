"""Tests for environment buffering (constant pH via fixed H+)."""

import math

import numpy as np

from ChemCompute import Compound, Enviroment, Reaction


def _aq(formula):
    return Compound(formula, phase_point_list=[{"phase": "aq", "temperature": 298}])


def _weak_acid_env(*, h_plus=1e-7, ha=0.1, a_minus=0.0, buffer=None):
    """HA <=> H+ + A- with buffered or free H+."""
    ha_c = _aq("HA")
    h_c = _aq("H+")
    a_c = _aq("A-")
    rxn = Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": ha_c, "rate_dependency": 1}],
        products=[
            {"stoichiometric_coefficient": 1, "compound": h_c, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": a_c, "rate_dependency": 1},
        ],
        reactants_concentration=[ha],
        products_concentration=[h_plus, a_minus],
        K=1e-5,
        kf=1.0,
        kb=1e5,
    )
    return Enviroment(
        rxn,
        concentrations={"H+": h_plus, "HA": ha, "A-": a_minus},
        buffer=buffer,
    )


class TestConstantPHEquilibrium:
    def test_buffered_h_plus_unchanged(self):
        env = _weak_acid_env(h_plus=1e-7, buffer=["H+"])
        result = env.equilibrium(return_details=True)
        assert np.isclose(result.concentrations_dict["H+"], 1e-7, rtol=1e-3)
        assert result.concentrations_dict["HA"] < 0.1
        assert result.concentrations_dict["A-"] > 0.0

    def test_unbuffered_h_plus_shifts(self):
        env = _weak_acid_env(h_plus=1e-7, buffer=None)
        result = env.equilibrium(return_details=True)
        assert result.concentrations_dict["H+"] > 1e-7

    def test_q_over_k_finite_with_buffered_h_plus(self):
        env = _weak_acid_env(h_plus=1e-7, buffer=["H+"])
        result = env.equilibrium(return_details=True)
        assert all(math.isfinite(q) for q in result.q_over_k)

    def test_explicit_buffer_target(self):
        env = _weak_acid_env(h_plus=1e-7, buffer={"H+": 1e-5})
        assert np.isclose(env.concentrations_dict["H+"], 1e-5)
        result = env.equilibrium(return_details=True)
        assert np.isclose(result.concentrations_dict["H+"], 1e-5, rtol=1e-3)


class TestConstantPHKinetics:
    def test_h_plus_flat_during_integration(self):
        env = _weak_acid_env(h_plus=1e-7, buffer=["H+"])
        traces = env.kinetics(time=1.0, accuracy=0.1)
        for snapshot in traces:
            h_index = env.compound_labels.index("H+")
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


class TestGenericBufferMechanism:
    def test_buffered_species_fixed_in_simple_equilibrium(self):
        a = _aq("A")
        b = _aq("B")
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
