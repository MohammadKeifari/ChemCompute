"""Tests for the Titration API."""

import math

import numpy as np

from ChemCompute import Compound, Enviroment, Reaction, Titration


def _water_env(h_plus=1e-7, oh_minus=1e-7, volume=0.1):
    h = Compound("H+", charge=1, phase_point_list=[{"temperature": 298, "phase": "aq"}])
    oh = Compound("OH-", charge=-1, phase_point_list=[{"temperature": 298, "phase": "aq"}])
    water = Compound("H2O", excess=True, phase_point_list=[{"temperature": 298, "phase": "l"}])
    rxn = Reaction(
        reactants=[
            {"stoichiometric_coefficient": 1, "compound": h, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": oh, "rate_dependency": 1},
        ],
        products=[{"stoichiometric_coefficient": 1, "compound": water, "rate_dependency": 0}],
        reactants_concentration=[h_plus, oh_minus],
        products_concentration=[0.0],
        K=1e-14,
    )
    return Enviroment(
        rxn,
        concentrations={"H+": h_plus, "OH-": oh_minus},
        volume=volume,
    )


class TestTitration:
    def test_strong_base_titration_pH_rises(self):
        sample = _water_env(h_plus=1e-3, oh_minus=1e-11, volume=0.1)
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
        sample = _water_env(volume=0.1)
        titrant = Enviroment.from_compounds({"OH-": 0.01, "Na+": 0.01}, volume=1.0)
        volumes = [0.0, 0.005, 0.01]
        result = Titration(sample, titrant, volumes=volumes).run(method="newton", tol=1e-8)

        assert result.titrant_volumes == volumes
        assert result.total_volumes[0] == sample.volume
        assert result.total_volumes[-1] > sample.volume

    def test_species_series(self):
        sample = _water_env(volume=0.1)
        titrant = Enviroment.from_compounds({"OH-": 0.01, "Na+": 0.01}, volume=1.0)
        result = Titration(sample, titrant, volume_min=0.0, volume_max=0.01, steps=5).run(
            method="newton",
            tol=1e-8,
        )

        h_series = result.species("H+")
        assert len(h_series) == 5
        assert np.all(np.diff(h_series) <= 0)

    def test_sample_not_mutated(self):
        sample = _water_env(h_plus=1e-3, volume=0.1)
        initial = sample.concentrations_dict.copy()
        titrant = Enviroment.from_compounds({"OH-": 0.1, "Na+": 0.1}, volume=1.0)
        Titration(sample, titrant, volume_min=0.0, volume_max=0.01, steps=5).run(
            method="newton",
            tol=1e-8,
        )
        assert sample.concentrations_dict == initial
        assert sample.volume == 0.1

    def test_temperature_mismatch_raises(self):
        sample = _water_env()
        titrant = Enviroment.from_compounds({"OH-": 0.1}, T=310, volume=1.0)
        try:
            Titration(sample, titrant)
            raise AssertionError("expected ValueError")
        except ValueError as exc:
            assert "temperature" in str(exc).lower()

    def test_plot_runs_without_display(self):
        sample = _water_env(volume=0.1)
        titrant = Enviroment.from_compounds({"OH-": 0.01, "Na+": 0.01}, volume=1.0)
        result = Titration(sample, titrant, volume_min=0.0, volume_max=0.01, steps=5).run(
            method="newton",
            tol=1e-8,
        )
        result.plot(species=["H+", "OH-"], plot="save", directory="titration_test_plot.png")
        result.plot_pH(plot="save", directory="titration_test_pH.png")

        import os

        assert os.path.isfile("titration_test_plot.png")
        assert os.path.isfile("titration_test_pH.png")
        os.remove("titration_test_plot.png")
        os.remove("titration_test_pH.png")
