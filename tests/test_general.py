import pytest
import numpy as np
from src.ChemCalc import Compound,Reaction,Enviroment


# -------------------------
#Compound Tests
# ------------------------

# ---------- Initialization & Unicode Handling ---------- #

def test_unicode_formula_generation_positive_charge():
    c = Compound("Na+1")
    # should contain superscript + and superscript ¹
    assert "\u207a" in c.unicode_formula
    assert "\u00b9" in c.unicode_formula


def test_unicode_formula_generation_negative_charge():
    c = Compound("Cl-1")
    # should contain superscript − and superscript ¹
    assert "\u207b" in c.unicode_formula
    assert "\u00b9" in c.unicode_formula


def test_unicode_formula_with_numbers_subscript():
    c = Compound("H2O")
    assert "\u2082" in c.unicode_formula  # ₂ for subscript 2
    assert "H" in c.unicode_formula
    assert "O" in c.unicode_formula


def test_scription_disabled_uses_plain_formula():
    c = Compound("CO2", scription=False)
    assert c.unicode_formula == "CO2"


# ---------- Phase Data Validation ---------- #

def test_valid_phase_point_list():
    phase_list = [{"phase": "s", "temperature": 0}, {"phase": "l", "temperature": 100}]
    c = Compound("H2O", phase_point_list=phase_list)
    assert len(c.phase_point_list) == 2


def test_invalid_phase_point_list_raises():
    bad_phase = [{"phase": "x", "temperature": 25}]
    with pytest.raises(ValueError):
        Compound("NaCl", phase_point_list=bad_phase)


# ---------- Phase Determination Logic ---------- #

def test_phase_exact_match_from_list():
    phase_list = [{"phase": "g", "temperature": 300}]
    c = Compound("CO2", phase_point_list=phase_list)
    assert c.phase(300) == "g"


def test_phase_solid_liquid_gas_ranges():
    c = Compound("H2O", mp=0, bp=100)
    assert c.phase(-5) == "s"  # below mp
    assert c.phase(50) == "l"  # between mp and bp
    assert c.phase(120) == "g"  # above bp


def test_phase_only_mp_defined():
    c = Compound("NaCl", mp=800, bp=None)
    assert c.phase(700) == "s"
    assert c.phase(900) == "l"


def test_phase_only_bp_defined():
    c = Compound("Ethanol", mp=None, bp=78)
    assert c.phase(20) == "l"
    assert c.phase(120) == "g"


def test_phase_no_mp_bp_returns_none():
    c = Compound("He")
    assert c.phase(300) is None


# ---------- Representation and Equality ---------- #

def test_str_returns_unicode_formula():
    c = Compound("CO2")
    assert str(c) == c.unicode_formula


def test_eq_compares_unicode_formulas():
    c1 = Compound("H2O")
    c2 = Compound("H2O")
    c3 = Compound("CO2")
    assert c1 == c2
    assert not (c1 == c3)


# ---------- Edge Cases ---------- #

def test_phase_point_list_and_mp_bp_combination():
    """Should prefer explicit temperature match in phase_point_list."""
    phase_list = [{"phase": "g", "temperature": 50}]
    c = Compound("NH3", phase_point_list=phase_list, mp=-78, bp=-33)
    assert c.phase(50) == "g"  # overrides mp/bp logic


def test_no_phase_points_returns_empty_list():
    c = Compound("O2")
    assert c.phase_point_list == []


# -------------------------
#Reaction Tests
# ------------------------

def test_reaction_initialization_basic():
    """Test direct initialization of a Reaction."""
    h2 = Compound("H2")
    o2 = Compound("O2")
    h2o = Compound("H2O")

    reactants = [{"stoichiometric_coefficient": 2, "compound": h2, "rate_dependency": 1}]
    products = [{"stoichiometric_coefficient": 2, "compound": h2o, "rate_dependency": 1}]

    r = Reaction(reactants, products, [1.0 , 1.0] , [0.0], K=2.0, kf=1.5, kb=0.5, T=300)

    assert len(r.reactants) == 1
    assert len(r.products) == 1
    assert r.K == 2.0
    assert r.kf == 1.5
    assert r.kb == 0.5
    assert isinstance(r.compounds, list)
    assert any(c["type"] == "reactant" for c in r.compounds)
    assert any(c["type"] == "product" for c in r.compounds)


def test_from_string_simple_syntax_parsing():
    """Check that a simple reaction string parses correctly."""
    reaction = Reaction.from_string_simple_syntax("2A.g + B.g2 > C.l-1", [1, 1, 0])
    assert len(reaction.reactants) == 2
    assert len(reaction.products) == 1
    assert isinstance(reaction.reactants[0]["compound"], Compound)
    assert reaction.reactants[0]["compound"].formula == "A"
    assert "⇌" in str(reaction)
    assert "(g)" in str(reaction)
    assert "(l)" in str(reaction)


def test_from_string_complex_syntax_parsing():
    """Check that a complex reaction string parses correctly."""
    reaction = Reaction.from_string_complex_syntax("2_H2_-1 & O2 > 2_H2O_-1", [1, 1, 0])
    assert isinstance(reaction, Reaction)
    assert len(reaction.reactants) == 2  # two reactants: H2 and O2
    assert len(reaction.products) == 1
    assert isinstance(reaction.products[0]["compound"], Compound)
    assert "H2" in reaction.products[0]["compound"].formula


def test_invalid_simple_syntax_raises():
    """Ensure invalid reaction strings raise ValueError."""
    with pytest.raises(ValueError):
        Reaction.from_string_simple_syntax("2H2 + + O2 > H2O")


def test_invalid_complex_syntax_raises():
    """Ensure invalid complex syntax raises ValueError."""
    with pytest.raises(ValueError):
        Reaction.from_string_complex_syntax("2H2 & > H2O")


def test_str_and_repr_output():
    """Check readable formatting."""
    reaction = Reaction.from_string_simple_syntax("A.g + B.g > C.g")
    text = repr(reaction)
    assert "⇌" in text
    assert "A" in text
    assert "(g)" in text


def test_addition_operator_combines_reactions():
    """Test combining two reactions with + operator."""
    r1 = Reaction.from_string_simple_syntax("A + B > C")
    r2 = Reaction.from_string_simple_syntax("2C > B + D")
    result = r1 + r2
    assert isinstance(result, Reaction)
    assert hasattr(result, "compounds")
    assert not("B" in repr(result))


def test_iadd_operator_aliases_add():
    """Test in-place addition behaves like +."""
    r1 = Reaction.from_string_simple_syntax("A + B > C")
    r2 = Reaction.from_string_simple_syntax("2C > B + D")
    combined = r1
    combined += r2
    assert isinstance(combined, Reaction)
    assert not("B" in repr(combined))


def test_iteration_over_compounds():
    """Ensure __iter__ yields compound dictionaries."""
    r = Reaction.from_string_simple_syntax("A.g + B.g > C.l")
    for item in r:
        assert "compound" in item
        assert "concentration" in item
        assert "type" in item


# -------------------------
#Enviroment Tests
# ------------------------
@pytest.fixture
def basic_env():
    """Create a small environment with a simple reversible reaction."""
    # A + 2B ⇌ C
    A = Compound("A")
    B = Compound("B")
    C = Compound("C")

    reactants = [
        {"stoichiometric_coefficient": 1, "compound": A, "rate_dependency": 1, "concentration": 1.0},
        {"stoichiometric_coefficient": 2, "compound": B, "rate_dependency": 2, "concentration": 2.0},
    ]
    products = [
        {"stoichiometric_coefficient": 1, "compound": C, "rate_dependency": 1, "concentration": 0.0}
    ]

    rxn = Reaction(reactants, products, [1.0 , 1.0], [0.0], K=10.0, kf=0.5, kb=0.1)
    env = Enviroment(rxn)
    return env


def test_check_if_reaction_accepts_only_reaction_objects(basic_env):
    """Check that _check_if_reaction validates correctly."""

    fake = "not_a_reaction"

    try:
        basic_env._check_if_reaction(fake)
    except ValueError:
        assert True
    else:
        assert False


def test_initialization_aggregates_compounds(basic_env):
    """Ensure the environment correctly aggregates compounds from its reactions."""
    compounds = [c.unicode_formula for c in basic_env.compounds]
    assert set(compounds) == {"A", "B", "C"}
    assert len(basic_env.compounds_concentration) == 3
    assert np.allclose(basic_env.concentrations_array, [1.0, 1.0, 0.0])


def test_add_reaction_via_method(basic_env):
    """Check adding a new reaction manually using add()."""
    D = Compound("D")
    E = Compound("E")
    new_rxn = Reaction(
        [{"stoichiometric_coefficient": 1, "compound": D, "rate_dependency": 1, "concentration": 1.5}],
        [{"stoichiometric_coefficient": 1, "compound": E, "rate_dependency": 1, "concentration": 0.2}],
        [1.0], [0.0], K=2.0, kf=0.3, kb=0.1,
    )
    initial_count = len(basic_env)
    basic_env.add(new_rxn)
    assert len(basic_env) == initial_count + 1
    assert D in basic_env.compounds
    assert E in basic_env.compounds


def test_add_reaction_with_iadd_operator(basic_env):
    """Check that += correctly adds a reaction."""
    D = Compound("D")
    E = Compound("E")
    new_rxn = Reaction(
        [{"stoichiometric_coefficient": 1, "compound": D, "rate_dependency": 1, "concentration": 0.5}],
        [{"stoichiometric_coefficient": 1, "compound": E, "rate_dependency": 1, "concentration": 0.5}],
        [1.0], [0.0], K=3.0, kf=0.4, kb=0.2,
    )
    initial = len(basic_env)
    basic_env += new_rxn
    assert len(basic_env) == initial + 1
    assert D in basic_env.compounds


def test_len_and_iter(basic_env):
    """Check __len__ and __iter__."""
    count = 0
    for rxn in basic_env:
        assert isinstance(rxn, Reaction)
        count += 1
    assert count == len(basic_env)


def test_reaction_by_index_mapping(basic_env):
    """Ensure that compound indices are correctly mapped in each reaction."""
    mapping = basic_env.reaction_by_index
    assert isinstance(mapping, list)
    assert all(isinstance(pair, list) for pair in mapping)
    reactant_idx, product_idx = mapping[0]
    assert len(reactant_idx) == 2
    assert len(product_idx) == 1


def test_stoichiometric_coefficient_array_shape_and_values(basic_env):
    """Verify correct stoichiometric matrix creation."""
    mat = basic_env.stoichiometric_coefficient_array
    assert mat.shape == (1, 3)
    # For A + 2B ⇌ C => [1, 2, -1]
    assert np.allclose(mat[0], [1, 2, -1])


def test_rate_constants_and_arrays(basic_env):
    """Ensure rate constants arrays are constructed correctly."""
    arr = basic_env.rate_constants_array
    assert arr.shape == (1, 2)
    assert np.allclose(arr[0], [0.5, 0.1])

    lst = basic_env.rate_constants
    assert isinstance(lst, list)
    assert lst[0] == [0.5, 0.1]


def test_rate_dependency_array_shape(basic_env):
    """Ensure rate dependency matrix is correct."""
    arr = basic_env.rate_dependency_array
    assert arr.shape == (1, 2, 3)
    # Reactants = [1,2,0], Products = [0,0,1]
    assert np.allclose(arr[0, 0], [1, 2, 0])
    assert np.allclose(arr[0, 1], [0, 0, 1])


def test_concentrations_get_and_set(basic_env):
    """Check concentration getter and setter behavior."""
    new_concs = [0.5, 1.0, 0.2]
    basic_env.concentrations = new_concs
    assert np.allclose(basic_env.concentrations_array, new_concs)

    # Wrong size should raise
    with pytest.raises(ValueError):
        basic_env.concentrations = [1.0, 2.0]


def test_unicode_formula_property(basic_env):
    """Verify compound Unicode formulas."""
    uforms = basic_env.compounds_unicode_formula
    assert isinstance(uforms, list)
    assert all(isinstance(u, str) for u in uforms)