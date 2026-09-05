"""Tests for expansion features: activity, thermo, buffer, scan, uvvis, bio."""

import math

import numpy as np

from ChemCompute import (
    ActivityModel,
    Compound,
    Enviroment,
    ParameterScan,
    Reaction,
    SpectrumSpec,
    buffer_diagnostics,
    competitive_inhibition,
    ionic_strength,
    single_substrate_mm,
    uvvis_spectrum,
)


class TestActivityModel:
    def test_davies_gamma_at_ionic_strength(self):
        env = Enviroment.__new__(Enviroment)
        na = Compound("Na+", charge=1, phase_point_list=[{"temperature": 298, "phase": "aq"}])
        cl = Compound("Cl-", charge=-1, phase_point_list=[{"temperature": 298, "phase": "aq"}])
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
        env.compounds = [
            Compound("Na+", charge=1, phase_point_list=[{"temperature": 298, "phase": "aq"}]),
            Compound("Cl-", charge=-1, phase_point_list=[{"temperature": 298, "phase": "aq"}]),
        ]
        env.charge_map = {}
        env._T = 298
        i = ionic_strength(env, np.array([0.2, 0.2]), 298)
        assert np.isclose(i, 0.2)


class TestVanHoffAndThermoFreeze:
    def test_full_vant_hoff_with_entropy(self):
        rxn = Reaction(
            reactants=[
                {
                    "stoichiometric_coefficient": 1,
                    "compound": Compound("A", phase_point_list=[{"temperature": 298, "phase": "aq"}]),
                    "rate_dependency": 1,
                }
            ],
            products=[
                {
                    "stoichiometric_coefficient": 1,
                    "compound": Compound("B", phase_point_list=[{"temperature": 298, "phase": "aq"}]),
                    "rate_dependency": 1,
                }
            ],
            reactants_concentration=[1.0],
            products_concentration=[0.0],
            K=1.0,
            enthalpy=5000.0,
            entropy=20.0,
            T=298,
        )
        rxn.T = 310
        manual = 1.0 * math.exp(
            -((5000.0 - 310.0 * 20.0) / (8.3145 * 310.0) - (5000.0 - 298.0 * 20.0) / (8.3145 * 298.0))
        )
        assert np.isclose(rxn.K, manual)

    def test_dh_only_vant_hoff_when_entropy_zero(self):
        rxn = Reaction(
            reactants=[
                {
                    "stoichiometric_coefficient": 1,
                    "compound": Compound("A", phase_point_list=[{"temperature": 298, "phase": "aq"}]),
                    "rate_dependency": 1,
                }
            ],
            products=[
                {
                    "stoichiometric_coefficient": 1,
                    "compound": Compound("B", phase_point_list=[{"temperature": 298, "phase": "aq"}]),
                    "rate_dependency": 1,
                }
            ],
            reactants_concentration=[1.0],
            products_concentration=[0.0],
            K=1.0,
            enthalpy=5000.0,
            entropy=0.0,
            T=298,
        )
        rxn.T = 310
        manual = 1.0 * math.exp(-5000.0 / 8.3145 * (1.0 / 310.0 - 1.0 / 298.0))
        assert np.isclose(rxn.K, manual)

    def test_adjust_thermodynamics_false_freezes_k(self):
        rxn = Reaction(
            reactants=[
                {
                    "stoichiometric_coefficient": 1,
                    "compound": Compound("A", phase_point_list=[{"temperature": 298, "phase": "aq"}]),
                    "rate_dependency": 1,
                }
            ],
            products=[
                {
                    "stoichiometric_coefficient": 1,
                    "compound": Compound("B", phase_point_list=[{"temperature": 298, "phase": "aq"}]),
                    "rate_dependency": 1,
                }
            ],
            reactants_concentration=[1.0],
            products_concentration=[0.0],
            K=2.5,
            kf=1.0,
            kb=0.4,
            enthalpy=10000.0,
            activation_energy_forward=5000.0,
            T=298,
        )
        env = Enviroment(rxn, adjust_thermodynamics=False)
        k0, kf0, kb0 = rxn.K, rxn.kf, rxn.kb
        env.T = 350
        assert rxn.K == k0
        assert rxn.kf == kf0
        assert rxn.kb == kb0


class TestBufferDiagnostics:
    def _ammonia_env(self):
        nh4 = Compound("NH4+", charge=1, phase_point_list=[{"temperature": 298, "phase": "aq"}])
        nh3 = Compound("NH3", phase_point_list=[{"temperature": 298, "phase": "aq"}])
        h = Compound("H+", charge=1, phase_point_list=[{"temperature": 298, "phase": "aq"}])
        oh = Compound("OH-", charge=-1, phase_point_list=[{"temperature": 298, "phase": "aq"}])
        water = Compound("H2O", excess=True)
        ka = 5.6e-10
        rxn1 = Reaction(
            reactants=[{"stoichiometric_coefficient": 1, "compound": nh4, "rate_dependency": 1}],
            products=[
                {"stoichiometric_coefficient": 1, "compound": h, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": nh3, "rate_dependency": 1},
            ],
            reactants_concentration=[0.05],
            products_concentration=[1e-9, 0.05],
            K=ka,
        )
        rxn2 = Reaction(
            reactants=[
                {"stoichiometric_coefficient": 1, "compound": h, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": oh, "rate_dependency": 1},
            ],
            products=[{"stoichiometric_coefficient": 1, "compound": water, "rate_dependency": 0}],
            reactants_concentration=[1e-9, 1e-7],
            products_concentration=[0.0],
            K=1e-14,
        )
        return Enviroment(rxn1, rxn2)

    def test_buffer_beta_positive(self):
        env = self._ammonia_env()
        result = env.equilibrium(method="newton", tol=1e-10, return_details=True)
        diag = buffer_diagnostics(env, result.concentrations)
        assert diag.beta > 0
        assert not math.isnan(diag.pH)


class TestParameterScan:
    def test_strong_acid_titration_pH_rises(self):
        h = Compound("H+", charge=1, phase_point_list=[{"temperature": 298, "phase": "aq"}])
        oh = Compound("OH-", charge=-1, phase_point_list=[{"temperature": 298, "phase": "aq"}])
        water = Compound("H2O", excess=True)
        kw = 1e-14
        rxn = Reaction(
            reactants=[
                {"stoichiometric_coefficient": 1, "compound": h, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": oh, "rate_dependency": 1},
            ],
            products=[{"stoichiometric_coefficient": 1, "compound": water, "rate_dependency": 0}],
            reactants_concentration=[1e-7, 1e-7],
            products_concentration=[0.0],
            K=kw,
        )
        env = Enviroment(rxn)
        scan = ParameterScan(
            base_env=env,
            axis="titrant_volume",
            titrant={"formula": "OH-", "concentration": 0.01, "volume_steps": np.linspace(0, 0.02, 20)},
            sample_volume=0.1,
        )
        curve = scan.run_equilibrium(method="newton", tol=1e-8)
        assert curve.pH[-1] > curve.pH[0]


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
        dye = Compound("D", phase_point_list=[{"temperature": 298, "phase": "aq"}])
        dummy = Compound("X", phase_point_list=[{"temperature": 298, "phase": "aq"}])
        rxn = Reaction(
            reactants=[
                {"stoichiometric_coefficient": 1, "compound": dye, "rate_dependency": 1},
            ],
            products=[
                {"stoichiometric_coefficient": 1, "compound": dummy, "rate_dependency": 1},
            ],
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
