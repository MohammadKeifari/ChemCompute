# Activity and buffering

## Ionic activity model

Correct Q/K at non-negligible ionic strength:

```python
from ChemCompute import ActivityModel

env = Enviroment(rxn, activity_model="davies")
env.charge_map = {"Ca+2": 2, "F-": -1}  # override auto-inferred charge
```

Supported models are listed in `VALID_ACTIVITY_MODELS`. Activity coefficients apply during `equilibrium()` when computing Q.

## Constant pH / buffering during solves

Hold pH fixed by setting `[H⁺]` and marking it buffered. Buffered species stay at their initial concentration while others react freely. Buffered species still appear in Q/K mass-action terms.

```python
env = Enviroment(
    weak_acid_rxn, water_autoionization,
    concentrations={"H+": 1e-7},
    buffer=["H+"],
)

env.apply_equilibrium(method="newton")
env.kinetics(time=10.0, accuracy=1e-3)  # [H+] stays 1e-7 M
```

For basic media: `concentrations={"OH-": 1e-2}, buffer=["OH-"]`.

Explicit targets: `buffer={"H+": 1e-7}`. Use `env.set_buffer(["H+"])` after changing concentrations.

!!! note "Buffering vs diagnostics"
    `buffer=[...]` **enforces** fixed concentration during the solve.
    `env.buffer_diagnostics()` only **reports** buffer capacity β(pH) after a calculation — it does not fix pH.

## Buffer diagnostics

After equilibrium:

```python
from ChemCompute import buffer_diagnostics

diag = env.buffer_diagnostics()
# diag.pH, diag.beta, diag.buffer_pairs, diag.hh_predictions
```

## Excess species (XS)

`XS(amount)` on reaction or environment concentration entries fixes amount during equilibrium/kinetics. Solids/liquids (`.s` / `.l`) are omitted from Q by phase regardless.

## Related

- [Equilibrium](../core/equilibrium.md)
- [Environment](../core/environment.md)
