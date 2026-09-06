# Env 16 equilibrium

HF / fluoride speciation in 0.06 M HCl with excess CaF₂(s) and H₂O(l). Excess solids and liquid water stay at activity 1; HCl is fully dissociated (`infinite_K`).

This network is also defined as `build_env16()` in `tests/manual/complex_environments.py`.

## Code

```python
from ChemCompute import Enviroment, Reaction, XS

env = Enviroment(
    Reaction.from_string("HF.aq & F-.aq > HF2-.aq", {}, K=0.1),
    Reaction.from_string("2_HF.aq > H2F2.aq", {}, K=0.5),
    Reaction.from_string("HF.aq > H+ & F-", {}, K=10 ** (-2.93)),
    Reaction.from_string(
        "H2O.l > H+ & OH-",
        {"H2O": XS(55.5), "H+": 1e-14 / 0.06},
        K=1e-14,
    ),
    Reaction.from_string("CaF2.s > Ca+2 & 2_F-", {"CaF2": XS(10.0)}, K=5e-9),
    Reaction.from_string("HCl.aq > H+ & Cl-", {"HCl": 0.06}, K=1.0, infinite_K=True),
)

result = env.equilibrium(method="newton", tol=1e-10, return_details=True)
print(result.concentrations_dict)
```

Or reuse the test builder:

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path("tests/manual")))
from complex_environments import build_env16

result = build_env16().equilibrium(method="newton", tol=1e-10, return_details=True)
```

## Figure

Bar chart of aqueous species from `result.concentrations_dict` (solids and water omitted):

![Env 16 aqueous equilibrium speciation](../images/env16_speciation.png)

## Related

- [Equilibrium](../core/equilibrium.md)
- [Reaction](../core/reaction.md)
- [Testing](../project/testing.md)
