"""
Manual validation for ChemCompute environments.

Run from the project root:

    python tests/manual/manual_validation.py

Add or edit environments in the section below (env1, env2, ...).
Each case runs equilibrium (printed to the terminal) and kinetics
(saves a plot under manual_test_output/kinetics/).
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ChemCompute import Compound, Enviroment, EquilibriumResult, Reaction

sys.path.insert(0, str(Path(__file__).resolve().parent))
from complex_environments import (
    COMPLEX_EQUILIBRIUM_KWARGS,
    ENV16_EXPECTED,
    ENV17_EXPECTED,
    ENV18_EXPECTED,
    ENV19_EXPECTED,
    ENV20_EXPECTED,
    build_env16,
    build_env17,
    build_env18,
    build_env19,
    build_env20,
)

OUTPUT_DIR = Path(__file__).resolve().parent.parent.parent / "manual_test_output"
KINETICS_DIR = OUTPUT_DIR / "kinetics"


@dataclass
class ManualCase:
    name: str
    description: str
    env: Enviroment
    expected_equilibrium: list[float]
    acceptable_rel_error: float
    equilibrium_kwargs: dict = field(default_factory=dict)
    run_kinetics: bool = True
    kinetic_kwargs: dict = field(default_factory=lambda: {"time": 5.0, "accuracy": 0.01})


# ---------------------------------------------------------------------------
# Environments (env1, env2, ...)
# ---------------------------------------------------------------------------

# env1: simple reversible reaction, batch gradient descent, log_quotient loss
env1 = Enviroment(
    Reaction.from_string_simple_syntax(
        "A.g > B.g",
        concentrations=[1.0, 0.0],
        K=2.0,
        kf=0.5,
        kb=0.25,
    ),
    T=298,
)

# env2: coupled network A ⇌ B ⇌ C, stochastic gradient descent
A2, B2, C2 = Compound("A"), Compound("B"), Compound("C")
env2 = Enviroment(
    Reaction(
        [{"stoichiometric_coefficient": 1, "compound": A2, "rate_dependency": 1}],
        [{"stoichiometric_coefficient": 1, "compound": B2, "rate_dependency": 1}],
        [1.0],
        [0.0],
        K=2.0,
        kf=0.5,
        kb=0.25,
    ),
    Reaction(
        [{"stoichiometric_coefficient": 1, "compound": B2, "rate_dependency": 1}],
        [{"stoichiometric_coefficient": 1, "compound": C2, "rate_dependency": 1}],
        [0.0],
        [0.0],
        K=1.5,
        kf=0.3,
        kb=0.2,
    ),
    T=298,
)

# env3: non-unity stoichiometry A + 2B ⇌ C
A3, B3, C3 = Compound("A"), Compound("B"), Compound("C")
env3 = Enviroment(
    Reaction(
        [
            {"stoichiometric_coefficient": 1, "compound": A3, "rate_dependency": 1},
            {"stoichiometric_coefficient": 2, "compound": B3, "rate_dependency": 2},
        ],
        [{"stoichiometric_coefficient": 1, "compound": C3, "rate_dependency": 1}],
        [1.0, 2.0],
        [0.0],
        K=10.0,
        kf=0.5,
        kb=0.05,
    ),
    T=298,
)

# env4: phase handling (gas reactant, liquid/solid products excluded from Q)
A4 = Compound("A", phase_point_list=[{"phase": "g", "temperature": 298}])
B4 = Compound("B", phase_point_list=[{"phase": "l", "temperature": 298}])
C4 = Compound("C", phase_point_list=[{"phase": "s", "temperature": 298}])
env4 = Enviroment(
    Reaction(
        [{"stoichiometric_coefficient": 1, "compound": A4, "rate_dependency": 1}],
        [
            {"stoichiometric_coefficient": 1, "compound": B4, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": C4, "rate_dependency": 1},
        ],
        [1.0],
        [0.0, 0.0],
        K=5.0,
        kf=0.3,
        kb=0.06,
    ),
    T=298,
)

# env5: temperature-dependent equilibrium (van't Hoff / Arrhenius)
env5 = Enviroment(
    Reaction.from_string_simple_syntax(
        "A > B",
        concentrations=[1.0, 0.0],
        K=2.0,
        kf=0.5,
        kb=0.25,
        enthalpy=-30000.0,
        entropy=-50.0,
        T=298,
    ),
    T=350,
)

# env6: Newton's method on a simple system
env6 = Enviroment(
    Reaction.from_string_simple_syntax(
        "A > B",
        concentrations=[2.0, 1.0],
        K=2.0,
        kf=0.5,
        kb=0.25,
    ),
    T=298,
)

# env7: quotient_error loss function
env7 = Enviroment(
    Reaction.from_string_simple_syntax(
        "A > B",
        concentrations=[1.0, 0.0],
        K=4.0,
        kf=0.8,
        kb=0.2,
    ),
    T=298,
)

# env8: log_huber loss function
env8 = Enviroment(
    Reaction.from_string_simple_syntax(
        "A > B",
        concentrations=[1.0, 0.0],
        K=0.5,
        kf=0.2,
        kb=0.4,
    ),
    T=298,
)

# env9: quotient_error_limit stopping criterion
env9 = Enviroment(
    Reaction.from_string_simple_syntax(
        "A > B",
        concentrations=[1.0, 0.0],
        K=2.0,
        kf=0.5,
        kb=0.25,
    ),
    T=298,
)

# env10: complex syntax, aqueous phase, multi-reaction system
env10 = Enviroment(
    Reaction.from_string_complex_syntax(
        "A.aq & 2_B.aq > C.aq",
        concentrations=[1.0, 1.0, 0.0],
        K=3.0,
        kf=0.4,
        kb=0.13333333333333333,
    ),
    Reaction.from_string_simple_syntax(
        "C.aq > D.aq",
        concentrations=[0.0, 0.0],
        K=1.2,
        kf=0.2,
        kb=0.16666666666666666,
    ),
    T=298,
)

# env11: weak acid dissociation with trace H+ (1e-8) and A- (1e-9)
env11 = Enviroment(
    Reaction(
        [{"stoichiometric_coefficient": 1, "compound": Compound("HA", phase_point_list=[{"phase": "aq", "temperature": 298}]), "rate_dependency": 1}],
        [
            {"stoichiometric_coefficient": 1, "compound": Compound("H+", phase_point_list=[{"phase": "aq", "temperature": 298}]), "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": Compound("A-", phase_point_list=[{"phase": "aq", "temperature": 298}]), "rate_dependency": 1},
        ],
        [0.1],
        [1e-8, 1e-9],
        K=1e-5,
        kf=1.0,
        kb=1e5,
    ),
    T=298,
)

# env12: acid-base neutralization H+(aq) + OH-(aq) -> H2O(l), trace ions
env12 = Enviroment(
    Reaction(
        [
            {"stoichiometric_coefficient": 1, "compound": Compound("H+", phase_point_list=[{"phase": "aq", "temperature": 298}]), "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": Compound("OH-", phase_point_list=[{"phase": "aq", "temperature": 298}]), "rate_dependency": 1},
        ],
        [{"stoichiometric_coefficient": 1, "compound": Compound("H2O", phase_point_list=[{"phase": "l", "temperature": 298}]), "rate_dependency": 1}],
        [1e-7, 1e-8],
        [55.0],
        K=1e14,
        kf=1e14,
        kb=1.0,
    ),
    T=298,
)

# env13: precipitation Ag+(aq) + Cl-(aq) -> AgCl(s)
env13 = Enviroment(
    Reaction(
        [
            {"stoichiometric_coefficient": 1, "compound": Compound("Ag+", phase_point_list=[{"phase": "aq", "temperature": 298}]), "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": Compound("Cl-", phase_point_list=[{"phase": "aq", "temperature": 298}]), "rate_dependency": 1},
        ],
        [{"stoichiometric_coefficient": 1, "compound": Compound("AgCl", phase_point_list=[{"phase": "s", "temperature": 298}]), "rate_dependency": 1}],
        [0.05, 0.05],
        [0.0],
        K=1e4,
        kf=1e4,
        kb=1.0,
    ),
    T=298,
)

# env14: dilute acid (1 mM) with ultra-trace ions (1e-9 / 1e-8)
env14 = Enviroment(
    Reaction(
        [{"stoichiometric_coefficient": 1, "compound": Compound("HA", phase_point_list=[{"phase": "aq", "temperature": 298}]), "rate_dependency": 1}],
        [
            {"stoichiometric_coefficient": 1, "compound": Compound("H+", phase_point_list=[{"phase": "aq", "temperature": 298}]), "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": Compound("A-", phase_point_list=[{"phase": "aq", "temperature": 298}]), "rate_dependency": 1},
        ],
        [1e-3],
        [1e-9, 1e-8],
        K=1e-4,
        kf=1.0,
        kb=1e4,
    ),
    T=298,
)

# env15: solid precipitation with pre-existing solid (Ag2CrO4)
env15 = Enviroment(
    Reaction(
        [
            {"stoichiometric_coefficient": 2, "compound": Compound("Ag+", phase_point_list=[{"phase": "aq", "temperature": 298}]), "rate_dependency": 2},
            {"stoichiometric_coefficient": 1, "compound": Compound("CrO4-2", phase_point_list=[{"phase": "aq", "temperature": 298}]), "rate_dependency": 1},
        ],
        [{"stoichiometric_coefficient": 1, "compound": Compound("Ag2CrO4", phase_point_list=[{"phase": "s", "temperature": 298}]), "rate_dependency": 1}],
        [0.2, 0.1],
        [0.5],
        K=1e6,
        kf=1e6,
        kb=1.0,
    ),
    T=298,
)


MANUAL_CASES: list[ManualCase] = [
    ManualCase(
        name="env1",
        description="Simple A <=> B (gas), BGD, log_quotient",
        env=env1,
        expected_equilibrium=[1.0 / 3.0, 2.0 / 3.0],
        acceptable_rel_error=0.02,
        equilibrium_kwargs={"method": "bgd", "loss": "log_quotient", "tol": 1e-8, "max_iter": 5000},
        kinetic_kwargs={"time": 8.0, "accuracy": 0.02},
    ),
    ManualCase(
        name="env2",
        description="Coupled A <=> B <=> C network, SGD",
        env=env2,
        expected_equilibrium=[1.0 / 6.0, 1.0 / 3.0, 0.5],
        acceptable_rel_error=0.03,
        equilibrium_kwargs={"method": "sgd", "loss": "log_quotient", "tol": 1e-8, "max_iter": 8000},
        kinetic_kwargs={"time": 10.0, "accuracy": 0.02},
    ),
    ManualCase(
        name="env3",
        description="Non-unity stoichiometry A + 2B <=> C",
        env=env3,
        expected_equilibrium=[0.26400109360177704, 0.5280021872035541, 0.735998906398223],
        acceptable_rel_error=0.02,
        equilibrium_kwargs={"method": "bgd", "loss": "log_quotient", "tol": 1e-8, "max_iter": 5000},
        kinetic_kwargs={"time": 12.0, "accuracy": 0.02},
    ),
    ManualCase(
        name="env4",
        description="Phase exclusion (g/l/s)",
        env=env4,
        expected_equilibrium=[0.2, 0.8, 0.8],
        acceptable_rel_error=0.02,
        equilibrium_kwargs={"method": "bgd", "loss": "log_quotient", "tol": 1e-8, "max_iter": 5000},
        kinetic_kwargs={"time": 6.0, "accuracy": 0.02},
    ),
    ManualCase(
        name="env5",
        description="Temperature-dependent K at 350 K",
        env=env5,
        expected_equilibrium=[0.7513342374511447, 0.24866576254885533],
        acceptable_rel_error=0.02,
        equilibrium_kwargs={"method": "bgd", "loss": "log_quotient", "tol": 1e-8, "max_iter": 5000},
        kinetic_kwargs={"time": 8.0, "accuracy": 0.02},
    ),
    ManualCase(
        name="env6",
        description="Newton's method, A <=> B with non-standard initial state",
        env=env6,
        expected_equilibrium=[1.0, 2.0],
        acceptable_rel_error=0.02,
        equilibrium_kwargs={"method": "newton", "loss": "log_quotient", "tol": 1e-10, "max_iter": 200},
        kinetic_kwargs={"time": 6.0, "accuracy": 0.01},
    ),
    ManualCase(
        name="env7",
        description="quotient_error loss",
        env=env7,
        expected_equilibrium=[0.2, 0.8],
        acceptable_rel_error=0.02,
        equilibrium_kwargs={"method": "bgd", "loss": "quotient_error", "tol": 1e-8, "max_iter": 5000},
        kinetic_kwargs={"time": 8.0, "accuracy": 0.02},
    ),
    ManualCase(
        name="env8",
        description="log_huber loss",
        env=env8,
        expected_equilibrium=[2.0 / 3.0, 1.0 / 3.0],
        acceptable_rel_error=0.03,
        equilibrium_kwargs={"method": "bgd", "loss": "log_huber", "huber_delta": 1.0, "tol": 1e-8, "max_iter": 5000},
        kinetic_kwargs={"time": 8.0, "accuracy": 0.02},
    ),
    ManualCase(
        name="env9",
        description="quotient_error_limit stopping (1%)",
        env=env9,
        expected_equilibrium=[1.0 / 3.0, 2.0 / 3.0],
        acceptable_rel_error=0.02,
        equilibrium_kwargs={
            "method": "bgd",
            "loss": "log_quotient",
            "quotient_error_limit": 0.01,
            "max_iter": 5000,
        },
        kinetic_kwargs={"time": 8.0, "accuracy": 0.02},
    ),
    ManualCase(
        name="env10",
        description="Complex syntax + aqueous multi-reaction network",
        env=env10,
        expected_equilibrium=[
            0.644535045447523,
            0.2890700908950461,
            0.16157497934203496,
            0.193889975210442,
        ],
        acceptable_rel_error=0.05,
        equilibrium_kwargs={"method": "newton", "loss": "log_quotient", "tol": 1e-10, "max_iter": 300},
        kinetic_kwargs={"time": 15.0, "accuracy": 0.02},
    ),
    ManualCase(
        name="env11",
        description="Weak acid HA(aq) with trace H+ (1e-8) and A- (1e-9)",
        env=env11,
        expected_equilibrium=[0.09900499297256835, 0.000995017027431657, 0.000995008027431657],
        acceptable_rel_error=0.02,
        equilibrium_kwargs={"method": "newton", "loss": "log_quotient", "tol": 1e-12, "max_iter": 200},
        kinetic_kwargs={"time": 8.0, "accuracy": 0.02},
    ),
    ManualCase(
        name="env12",
        description="Acid-base neutralization H+(aq) + OH-(aq) -> H2O(l), trace ions",
        env=env12,
        expected_equilibrium=[1.5465856099730644e-07, 6.465856099730644e-08, 54.999999945341436],
        acceptable_rel_error=0.05,
        equilibrium_kwargs={"method": "newton", "loss": "log_quotient", "tol": 1e-12, "max_iter": 200},
        kinetic_kwargs={"time": 1e-6, "accuracy": 1e-10},
    ),
    ManualCase(
        name="env13",
        description="Precipitation Ag+(aq) + Cl-(aq) -> AgCl(s)",
        env=env13,
        expected_equilibrium=[0.009999999999999988, 0.009999999999999988, 0.040000000000000015],
        acceptable_rel_error=0.02,
        equilibrium_kwargs={"method": "newton", "loss": "log_quotient", "tol": 1e-12, "max_iter": 200},
        kinetic_kwargs={"time": 6.0, "accuracy": 0.02},
    ),
    ManualCase(
        name="env14",
        description="Dilute acid (1 mM) with ultra-trace ions (1e-9 / 1e-8)",
        env=env14,
        expected_equilibrium=[0.0007298484291421942, 0.00027015257085780586, 0.00027016157085780584],
        acceptable_rel_error=0.02,
        equilibrium_kwargs={"method": "newton", "loss": "log_quotient", "tol": 1e-12, "max_iter": 200},
        kinetic_kwargs={"time": 8.0, "accuracy": 0.02},
    ),
    ManualCase(
        name="env15",
        description="Solid Ag2CrO4(s) precipitation with pre-existing solid",
        env=env15,
        expected_equilibrium=[0.012599210498948743, 0.0062996052494743715, 0.5937003947505256],
        acceptable_rel_error=0.02,
        equilibrium_kwargs={"method": "newton", "loss": "log_quotient", "tol": 1e-12, "max_iter": 200},
        kinetic_kwargs={"time": 6.0, "accuracy": 0.02},
    ),
    ManualCase(
        name="env16",
        description="HF fluoride speciation in 0.06 M HCl with excess CaF2(s) and H2O(l)",
        env=build_env16(),
        expected_equilibrium=ENV16_EXPECTED,
        acceptable_rel_error=0.05,
        equilibrium_kwargs=COMPLEX_EQUILIBRIUM_KWARGS["env16"],
        run_kinetics=False,
    ),
    ManualCase(
        name="env17",
        description="Carbonate buffering with excess CaCO3(s) and 0.01 M HCl (infinite K)",
        env=build_env17(),
        expected_equilibrium=ENV17_EXPECTED,
        acceptable_rel_error=0.05,
        equilibrium_kwargs=COMPLEX_EQUILIBRIUM_KWARGS["env17"],
        run_kinetics=False,
    ),
    ManualCase(
        name="env18",
        description="NH4+ buffer, AgCl(s) precipitation, HCl dissociation (infinite K)",
        env=build_env18(),
        expected_equilibrium=ENV18_EXPECTED,
        acceptable_rel_error=0.05,
        equilibrium_kwargs=COMPLEX_EQUILIBRIUM_KWARGS["env18"],
        run_kinetics=False,
    ),
    ManualCase(
        name="env19",
        description="Phosphate buffer with CaHPO4(s) precipitation and excess solid",
        env=build_env19(),
        expected_equilibrium=ENV19_EXPECTED,
        acceptable_rel_error=0.05,
        equilibrium_kwargs=COMPLEX_EQUILIBRIUM_KWARGS["env19"],
        run_kinetics=False,
    ),
    ManualCase(
        name="env20",
        description="Cu-EDTA complexation (infinite K), NaCl dissociation, water autoionization",
        env=build_env20(),
        expected_equilibrium=ENV20_EXPECTED,
        acceptable_rel_error=0.05,
        equilibrium_kwargs=COMPLEX_EQUILIBRIUM_KWARGS["env20"],
        run_kinetics=False,
    ),
]


def _relative_errors(computed: list[float], expected: list[float], floor: float = 1e-12) -> list[float]:
    errors = []
    for calc, exp in zip(computed, expected):
        denom = max(abs(exp), floor)
        errors.append(abs(calc - exp) / denom)
    return errors


def _within_tolerance(computed: list[float], expected: list[float], acceptable_rel_error: float) -> bool:
    errors = _relative_errors(computed, expected)
    return max(errors) <= acceptable_rel_error


def _format_concentration(value: float) -> str:
    if value == 0.0:
        return "0.0"
    if abs(value) >= 1e-3 or abs(value) < 1e-6:
        return f"{value:.6g}"
    return f"{value:.6e}"


def _format_quotient_error(value: float) -> str:
    if value == 0.0:
        return "0.0"
    if abs(value) >= 1e-4:
        return f"{value:.6f}"
    return f"{value:.6e}"


def _print_equilibrium_report(case: ManualCase, result: EquilibriumResult) -> bool:
    env = case.env
    computed = result.concentrations
    errors = _relative_errors(computed, case.expected_equilibrium)
    passed = _within_tolerance(computed, case.expected_equilibrium, case.acceptable_rel_error)

    print(f"\n{'=' * 72}")
    print(f"{case.name}: {case.description}")
    print(f"{'=' * 72}")
    print(f"Temperature (K): {env.T}")
    print(f"Compounds: {result.compounds}")
    print(f"Initial concentrations: {env.concentrations}")
    for idx, rxn in enumerate(env.reactions, start=1):
        print(f"Reaction {idx} K: {rxn.K:.6g}  kf: {rxn.kf:.6g}  kb: {rxn.kb:.6g}")
    print(f"Equilibrium settings: {case.equilibrium_kwargs}")
    print(f"Expected equilibrium:   {[_format_concentration(v) for v in case.expected_equilibrium]}")
    print(f"Computed equilibrium:   {[_format_concentration(v) for v in computed]}")
    print(f"Relative errors:        {[round(e, 6) for e in errors]}")
    print(f"Reaction extents (x):   {[round(v, 6) for v in result.reaction_extents]}")
    print(f"|Q/K - 1| per rxn:    {[_format_quotient_error(v) for v in result.reaction_quotient_error]}")
    print(f"Max |Q/K - 1|:        {_format_quotient_error(result.max_reaction_quotient_error)}")
    print(f"Q/K ratios:             {[_format_concentration(v) for v in result.reaction_quotient_ratio]}")
    if "H+" in result.concentrations_dict and "OH-" in result.concentrations_dict:
        kw = result.concentrations_dict["H+"] * result.concentrations_dict["OH-"]
        print(f"Kw check [H+][OH-]:   {kw:.6e}  (target 1e-14)")
    print(f"Criterion ({result.criterion_type}): {result.criterion_value:.6g} / {result.criterion_limit:.6g} -> {'MET' if result.criterion_met else 'NOT MET'}")
    print(f"Stop reason:            {result.stop_reason} ({result.iterations} iterations)")
    print(f"Acceptable max error:   {case.acceptable_rel_error * 100:.1f}%")
    print(f"Result:                 {'PASS' if passed else 'FAIL'}")
    return passed


def _run_kinetics(case: ManualCase) -> Path:
    KINETICS_DIR.mkdir(parents=True, exist_ok=True)
    plot_path = KINETICS_DIR / f"{case.name}.png"
    case.env.kinetics(
        plot="save",
        directory=str(plot_path),
        **case.kinetic_kwargs,
    )
    return plot_path


def run_manual_validation() -> int:
    print("ChemCompute manual validation")
    print(f"Kinetic plots output directory: {KINETICS_DIR}")

    passed_count = 0
    results: list[tuple[str, bool]] = []

    for case in MANUAL_CASES:
        equilibrium_kwargs = {**case.equilibrium_kwargs, "return_details": True}
        result = case.env.equilibrium(**equilibrium_kwargs)
        passed = _print_equilibrium_report(case, result)
        results.append((case.name, passed))
        if passed:
            passed_count += 1

        if case.run_kinetics:
            plot_path = _run_kinetics(case)
            print(f"Kinetics plot saved: {plot_path}")

    total = len(MANUAL_CASES)
    print(f"\n{'=' * 72}")
    print("SUMMARY")
    print(f"{'=' * 72}")
    for name, passed in results:
        print(f"  {name}: {'PASS' if passed else 'FAIL'}")
    print(f"\nEnvironments within acceptable error rate: {passed_count}/{total}")
    print(f"{'=' * 72}\n")

    return 0 if passed_count == total else 1


if __name__ == "__main__":
    raise SystemExit(run_manual_validation())
