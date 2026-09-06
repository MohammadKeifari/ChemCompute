# Equilibrium

`env.equilibrium()` finds concentrations that satisfy all finite-K reactions simultaneously by minimizing a loss on ln(Q/K) (or related metrics) subject to non-negative concentrations.

## Basic usage

```python
equilibrium = env.equilibrium(method="newton", tol=1e-10)

result = env.equilibrium(method="newton", tol=1e-10, return_details=True)
print(result.concentrations_dict)
print(result.criterion_met)
print(result.stop_reason)
```

Apply the solution onto the environment:

```python
result = env.apply_equilibrium(method="newton", tol=1e-10)
# env.concentrations updated; details in env.last_equilibrium_result
```

Default method is **`newton`** (best for stiff speciation networks).

## Methods and stopping criteria

| Parameter | Options | Role |
|-----------|---------|------|
| `method` | `"bgd"`, `"sgd"`, `"newton"` | Optimizer |
| `loss` | `"log_quotient"`, `"quotient_error"`, `"log_huber"` | Objective |
| `tol` | float | Residual tolerance (when `quotient_error_limit` is not set) |
| `quotient_error_limit` | float, e.g. `0.01` | Stop when every reaction has \|Q/K − 1\| ≤ limit |
| `max_iter`, `learning_rate` | — | Iteration control |
| `min_concentration` | float | Floor for log safety |

## Diagnostics

With `return_details=True` or via `env.last_equilibrium_result`:

- `reaction_quotient_error` — per-reaction \|Q/K − 1\|
- `reaction_quotient_ratio` — per-reaction Q/K
- `reaction_extents` — stoichiometric extent vector
- `criterion_met`, `criterion_type`, `criterion_value`, `criterion_limit`

## Phases, Kw, and Ksp

Solid (`s`) and liquid (`l`) species are **omitted from Q** with **activity = 1**:

- **Kw** = [H⁺][OH⁻] = K — H₂O does not appear in Q or K
- **Ksp** = [Ca²⁺][F⁻]² — the pure solid does not appear in Q or K

Only `aq` and `g` species participate in the mass-action product Q.

```python
from ChemCompute import XS

kw = Reaction.from_string(
    "H2O.l > H+ & OH-",
    concentrations={"H2O": XS(55.5), "H+": 1e-7, "OH-": 1e-7},
    K=1e-14,
)
ksp = Reaction.from_string(
    "CaF2.s > Ca+2 & 2_F-",
    concentrations={"CaF2": XS(10.0)},
    K=5e-9,
)
env = Enviroment(kw, ksp)
```

`XS(amount)` fixes concentration during the solve without multiplying K by bulk molarity.

## Activity corrections

Pass `activity_model="davies"` (or other supported models) to the environment constructor. See [Activity and buffering](../guide/activity-and-buffering.md).

## Redox coupling

Half-reactions update linked reaction K values from Nernst electrode potential when `env.electrode_Eh` is set or when multiple half-reactions share a solved Eh.

## Reference environments

Complex test cases env16–env20 live in `tests/manual/complex_environments.py`. Regenerate reference concentrations:

```bash
python tests/manual/generate_expected.py
```

See [Env 16 example](../examples/env16-equilibrium.md).

## Related

- [Reaction](reaction.md)
- [Environment](environment.md)
- [Pourbaix guide](../guide/pourbaix.md) — equilibrium at every grid point (slow)
