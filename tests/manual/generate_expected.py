"""Generate reference equilibrium concentrations for env16–env20.

Run from project root:

    python tests/manual/generate_expected.py

Uses a high-accuracy Newton solve (strict tol, low min_concentration).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from complex_environments import (
    REFERENCE_EQUILIBRIUM_KWARGS,
    build_env16,
    build_env17,
    build_env18,
    build_env19,
    build_env20,
)

BUILDERS = {
    "env16": build_env16,
    "env17": build_env17,
    "env18": build_env18,
    "env19": build_env19,
    "env20": build_env20,
}


def main() -> int:
    for name, builder in BUILDERS.items():
        env = builder()
        result = env.equilibrium(**REFERENCE_EQUILIBRIUM_KWARGS, return_details=True)
        labels = [c.formula for c in env.compounds]
        print(f"# {name}: compounds = {labels!r}")
        print(f"{name.upper()}_EXPECTED = {result.concentrations!r}")
        print(
            f"# {name}: stop={result.stop_reason}, iter={result.iterations}, "
            f"criterion_met={result.criterion_met}, "
            f"max_q_err={result.max_reaction_quotient_error:.6g}, "
            f"residual={result.criterion_value:.6g}"
        )
        for label, value in zip(labels, result.concentrations):
            print(f"#   {label:8s} {value:.12e}")
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
