import pytest
from ChemCompute import Enviroment, Reaction, uvvis_spectrum
from ChemCompute.compounds import (
    all_compounds,
    fescn,
    fe3,
    get,
    h2o,
    h_plus,
    mno4,
    oh_minus,
    scn_minus,
    so4_2minus,
    water,
)


def test_water_alias_and_kelvin_phase_points():
    assert water is h2o
    assert water() is h2o()
    assert water().formula == "H2O"
    assert water().mp == 273.15
    assert water().bp == 373.15
    assert water().phase(298.0) == "l"
    assert water().phase(250.0) == "s"
    assert water().phase(400.0) == "g"


def test_token_interpolation_keeps_library_object():
    rxn = Reaction.from_string(
        f"{water().token} > {h_plus().token} & {oh_minus().token}",
        concentrations=[1.0, 1e-7, 1e-7],
        K=1e-14,
    )
    assert rxn.reactants[0]["compound"] is water()
    assert rxn.products[0]["compound"] is h_plus()
    assert rxn.products[1]["compound"] is oh_minus()


def test_ion_charge_and_formula_lookup():
    assert so4_2minus().charge == -2
    assert h_plus().charge == 1
    assert get("H2O") is water()
    assert get("SO4-2") is so4_2minus()
    assert get("MnO4-") is mno4()


def test_library_formulas_are_unique():
    by_formula = {}
    for compound in all_compounds():
        by_formula.setdefault(compound.formula, []).append(compound)
    clashes = {formula: items for formula, items in by_formula.items() if len(items) > 1}
    assert clashes == {}


def test_permanganate_envelope_connects_through_the_visible_band():
    spec = mno4().spectrum
    assert spec.extrapolate == "flat"
    assert len(spec.points) > 3
    peak = spec.epsilon(525 * 1e-9)
    assert peak == pytest.approx(2450.0 * 100.0)
    assert 0.0 < spec.epsilon(535 * 1e-9) < peak
    assert spec.epsilon(700 * 1e-9) == 0.0


def test_fescn_absorbance_at_tabulated_peak():
    rxn = Reaction.from_string(
        f"{fe3().token} & {scn_minus().token} > {fescn().token}",
        concentrations=[0.0, 0.0, 1.0e-4],
        K=1.0,
    )
    env = Enviroment(rxn)
    absorbance = uvvis_spectrum(env, wavelengths=[447 * 1e-9], path_length=0.01)[0]
    assert absorbance == pytest.approx(4700.0 * 1.0e-4)
    off_peak = uvvis_spectrum(env, wavelengths=[470 * 1e-9], path_length=0.01)[0]
    assert 0.0 < off_peak < absorbance


def test_library_spectra_are_connected_envelopes():
    with_spec = [c for c in all_compounds() if c.spectrum is not None]
    assert with_spec
    for compound in with_spec:
        assert compound.spectrum.extrapolate == "flat"
        assert len(compound.spectrum.points) >= 5
        wavelengths = [point[0] for point in compound.spectrum.points]
        assert wavelengths[0] < wavelengths[-1]
        assert compound.spectrum.points[0][1] == 0.0
        assert compound.spectrum.points[-1][1] == 0.0
