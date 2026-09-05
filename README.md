# ChemCompute

ChemCompute models multi-reaction chemical systems in Python. You define compounds and reactions once, wrap them in an `Enviroment`, and then run either:

- **Equilibrium** — solve for concentrations where each reaction satisfies its mass-action expression (Q/K)
- **Kinetics** — integrate concentrations forward in time from rate laws

Both paths share the same reaction network, stoichiometry, and concentration state.

## Installation

```bash
pip install chemcompute
```

From source:

```bash
git clone <repository-url>
cd ChemCompute
pip install -e .
```

Requires Python 3.7+, NumPy, and Matplotlib (for plots).

```python
from ChemCompute import Compound, Reaction, Enviroment, EquilibriumResult
```

---

## Build a reaction system

Every calculation starts with compounds, reactions, and an environment.

```python
from ChemCompute import Compound, Reaction, Enviroment

# Simple syntax: A <=> B
rxn = Reaction.from_string_simple_syntax(
    "A > B",
    concentrations=[1.0, 0.0],
    K=2.0,
    kf=0.5,
    kb=0.25,
)

env = Enviroment(rxn, T=298)  # Kelvin
```

Multiple reactions share compounds automatically:

```python
rxn1 = Reaction.from_string_simple_syntax("A > B", [1.0, 0.0], K=2.0, kf=0.5, kb=0.25)
rxn2 = Reaction.from_string_simple_syntax("B > C", [0.0, 0.0], K=1.5, kf=0.3, kb=0.2)

env = Enviroment(rxn1, rxn2, T=298)
env.concentrations = [1.0, 0.0, 0.0]  # [A, B, C]
```

Useful accessors:

| Property / method | Purpose |
|-------------------|---------|
| `env.compounds` | Ordered list of unique species |
| `env.concentrations` | Current molar concentrations |
| `env.concentrations_dict` | `{formula: concentration}` mapping |
| `env.compound_labels` | Formula strings aligned with concentration vectors |

Reactions can also be built with `Reaction.from_string_complex_syntax()` or explicit reactant/product lists. See `tests/manual/manual_validation.py` for varied examples (phases, buffers, precipitation, coupled networks).

---

## Equilibrium calculation

`env.equilibrium()` finds concentrations that satisfy all finite-K reactions simultaneously. It minimizes a loss on ln(Q/K) (or related metrics) subject to non-negative concentrations.

### Basic usage

```python
# Returns a concentration list (same order as env.compounds)
equilibrium = env.equilibrium(method="newton", tol=1e-10)

# Or return full diagnostics
result = env.equilibrium(method="newton", tol=1e-10, return_details=True)
print(result.concentrations_dict)
print(result.criterion_met)   # Did the solver meet its stopping criterion?
print(result.stop_reason)     # e.g. "residual_tol", "quotient_error_limit"
```

Apply the solution back onto the environment:

```python
result = env.apply_equilibrium(method="newton", tol=1e-10)
# env.concentrations is now updated; details in env.last_equilibrium_result
```

### Methods and stopping criteria

| Parameter | Options | Role |
|-----------|---------|------|
| `method` | `"bgd"`, `"sgd"`, `"newton"` | Optimizer (Newton is usually best for stiff/speciation systems) |
| `loss` | `"log_quotient"`, `"quotient_error"`, `"log_huber"` | What to minimize |
| `tol` | float | Residual tolerance (default when `quotient_error_limit` is not set) |
| `quotient_error_limit` | float, e.g. `0.01` | Stop when every reaction has \|Q/K − 1\| ≤ limit (overrides `tol`) |
| `max_iter`, `learning_rate` | — | Iteration control (defaults depend on method) |
| `min_concentration` | float | Floor used only for log safety on zero/negative values |

**Diagnostics** (with `return_details=True` or via `env.last_equilibrium_result`):

- `reaction_quotient_error` — per-reaction \|Q/K − 1\|
- `reaction_quotient_ratio` — per-reaction Q/K (1.0 means equilibrium for that reaction)
- `reaction_extents` — stoichiometric extent vector (informational)
- `criterion_met`, `criterion_type`, `criterion_value`, `criterion_limit`

### Phases, Kw, and Ksp

Solid (`s`) and liquid (`l`) species are **omitted from Q** with **activity = 1**, matching standard thermodynamic practice:

- Water autoionization: **Kw = [H⁺][OH⁻] = K** — H₂O does not appear in Q or K
- Solid solubility: **Ksp = [Ca²⁺][F⁻]²** — the pure solid does not appear in Q or K

Only `aq` and `g` species participate in the mass-action product Q.

```python
h2o = Compound("H2O", phase_point_list=[{"phase": "l", "temperature": 298}], excess=True)
caf2 = Compound("CaF2", phase_point_list=[{"phase": "s", "temperature": 298}], excess=True)
```

`excess=True` keeps that species' concentration **fixed** during the solve (large reservoir of solid or solvent). It does not multiply K by the bulk molarity.

### Advanced reaction options

```python
# Strong acid fully dissociated — treated as irreversible (excluded from Q/K residual)
Reaction(..., infinite_K=True)

# Direct initialization with explicit Compound objects and non-unity stoichiometry
Reaction(reactants, products, reactants_concentration, products_concentration, K=..., kf=..., kb=...)
```

Complex reference environments (env16–env20) live in `tests/manual/complex_environments.py`. Regenerate their reference concentrations with:

```bash
python tests/manual/generate_expected.py
```

---

## Kinetic simulation

`env.kinetics()` integrates the mass-action rate laws forward in time using the current concentrations and rate constants `kf` / `kb`.

### Basic usage

```python
results = env.kinetics(
    time=10.0,           # total simulation time
    accuracy=1e-3,       # time step
    checkpoint_time=[1.0, 5.0, 10.0],
)

# results is a list of concentration snapshots at checkpoints (and final time)
final = results[-1]
```

### Plotting

```python
env.kinetics(
    time=10.0,
    plot="interactive",  # or "save" or False
    directory="./plot.png",
    colors=["#26547c", "#ef476f"],  # one color per compound; optional
)
```

Kinetics uses the same `env.concentrations` as the starting point. Run equilibrium first if you want to integrate from an equilibrated state:

```python
env.apply_equilibrium(method="newton", tol=1e-10)
env.kinetics(time=5.0, plot="save", directory="approach.png")
```

---

## Temperature-dependent K and k

Pass enthalpy, entropy, and activation energies when defining a reaction. Changing `rxn.T` or `env.T` updates **K** (van't Hoff) and **kf** / **kb** (Arrhenius):

```python
rxn = Reaction.from_string_simple_syntax(
    "A > B",
    K=2.0, kf=0.5, kb=0.25,
    enthalpy=-50000,
    activation_energy_forward=50000,
    activation_energy_backward=100000,
    T=298,
)
rxn.T = 350  # K, kf, kb recalculate automatically
```

---

## Testing

```bash
pytest tests/
python tests/manual/manual_validation.py   # 20 named equilibrium/kinetics cases
```

Kinetic plots from manual validation are written to `manual_test_output/kinetics/`.

---

## Project layout

```
src/ChemCompute/
  _general.py      Compound, Reaction, Enviroment
  _equilibrium.py  Equilibrium solver and EquilibriumResult
  _kinetics.py     Time integration and plotting
tests/
  test_environment_calculators.py
  manual/          Reference environments and validation scripts
docs/index.md      Extended documentation
```

---

## License

MIT — see [LICENSE](LICENSE).

## Author

Mohammad Keifari — [mohammadkeifari2007@gmail.com](mailto:mohammadkeifari2007@gmail.com)
