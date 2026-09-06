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

## Environment construction and mixing

### Three ways to build an environment

```python
# 1) Classic — reactions only (concentrations summed from reactions)
env = Enviroment(rxn1, rxn2, T=298, volume=1.0)

# 2) Compounds only — no reaction network
env = Enviroment.from_compounds({"Na+": 0.1, "Cl-": 0.1}, volume=1.0)

# 3) Reactions + concentration overrides (dict wins over reaction values)
env = Enviroment(
    rxn1, rxn2,
    concentrations={"H+": 0.06, "Ca+2": 0.01},
    volume=0.1,
)
```

`volume` is in litres (default `1.0`) and is used when mixing environments.

### Standalone reaction solvers

```python
rxn = Reaction.from_string_simple_syntax("A > B", [1.0, 0.0], K=2.0, kf=0.5, kb=0.25)
final = rxn.equilibrium(method="newton", tol=1e-10)
traces = rxn.kinetics(time=5.0, accuracy=1e-3)
```

### Mix environments (volume-weighted)

```python
envA = Enviroment(rxn_a, concentrations={"A": 1.0}, volume=1.0)
envB = Enviroment.from_compounds({"B": 0.1}, volume=1.0)

envC = envA + envB                         # volume = 2.0
envD = 0.5 * envA + 4 * envB               # effective volume = 4.5 L
envE = Enviroment.combine((0.5, envA), (4, envB))  # same as envD

envF = envC.add_compounds({"A": 1.0}, volume=1.0)
envG = envC.add_compounds({"A": 1.0}, volume=1.0, coefficient=4.0)
```

Mixing rule: `effective_volume = coeff × volume`, total moles per species are conserved, final concentration = moles / total effective volume.

### Constant pH / buffering

Hold pH fixed during equilibrium or kinetics by setting `[H⁺]` and marking it as buffered. ChemCompute keeps buffered species at their initial concentration (after `concentrations=` overrides) while other species react freely. Buffered species still appear in Q/K mass-action terms.

```python
env = Enviroment(
    weak_acid_rxn, water_autoionization,
    concentrations={"H+": 1e-7},   # pH 7
    buffer=["H+"],                 # hold [H+] constant during the solve
)

env.apply_equilibrium(method="newton")
env.kinetics(time=10.0, accuracy=1e-3)  # [H+] stays 1e-7 M
```

For basic media, fix `[OH⁻]` instead: `concentrations={"OH-": 1e-2}, buffer=["OH-"]`.

Explicit targets are optional: `buffer={"H+": 1e-7}`. Use `env.set_buffer(["H+"])` to re-snapshot after changing concentrations.

This **enforces** constant concentration during the solve. `env.buffer_diagnostics()` only **reports** buffer capacity β(pH) after a calculation—it does not fix pH. `Compound.excess=True` fixes solids/liquids and omits them from Q; buffered H⁺ stays in Q at its fixed value.

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

When `entropy` is set, full van't Hoff uses ΔG = ΔH − TΔS. Disable automatic updates:

```python
env = Enviroment(rxn, adjust_thermodynamics=False)
env.T = 350  # K, kf, kb unchanged
```

---

## Ionic activity model

Correct Q/K at non-negligible ionic strength with Debye–Hückel, Davies, or Pitzer-lite:

```python
from ChemCompute import ActivityModel

env = Enviroment(rxn, activity_model="davies")
env.charge_map = {"Ca2+": 2, "F-": -1}  # overrides Compound.charge
# or: Compound("Na+", charge=1, ...)
```

Activity coefficients are applied in `equilibrium()` when computing Q.

---

## Buffer diagnostics

After equilibrium, compute buffer capacity β(pH) and Henderson–Hasselbalch checks:

```python
from ChemCompute import buffer_diagnostics

diag = env.buffer_diagnostics()  # uses current concentrations
# diag.pH, diag.beta, diag.buffer_pairs, diag.hh_predictions
```

---

## Titration

Mix a **sample** environment with a **titrant** environment over a range of titrant volumes, equilibrate at each step, and collect concentrations vs volume added.

```python
from ChemCompute import Titration, Enviroment

sample = Enviroment(...)  # analyte; set sample.volume (e.g. 0.1 L)
titrant = Enviroment.from_compounds({"OH-": 0.1, "Na+": 0.1}, volume=1.0)

curve = Titration(
    sample,
    titrant,
    volume_min=0.0,
    volume_max=0.05,
    steps=100,
).run(method="newton", tol=1e-10)

# Access data
curve.titrant_volumes          # L added at each step
curve.pH
curve.matrix()                 # shape (n_steps, n_compounds)
curve.species("H+")            # one species vs volume
curve.speciation               # list of {formula: c} per step

# Plot
curve.plot(species=["H+", "OH-"], plot="save", directory="titration.png")
curve.plot_pH(plot="interactive")
```

Pass an explicit volume list with ``volumes=[0, 0.001, 0.002, ...]`` instead of ``volume_min``/``volume_max``/``steps``.

The sample and titrant are not mutated. Mixing uses the same volume-weighted rules as ``envA + envB``.

---

## UV–Vis (Beer–Lambert)

User-supplied molar absorptivity points with piecewise-linear interpolation:

```python
from ChemCompute import SpectrumSpec, uvvis_spectrum

env.set_spectrum("InH", SpectrumSpec(
    points=[(400e-9, 12000), (450e-9, 25000), (500e-9, 8000)],
    extrapolate="flat",  # or "none" for zero epsilon off-tabulated wavelengths
))
A = uvvis_spectrum(env, wavelengths=[450e-9, 500e-9], path_length=0.01)
```

---

## Enzyme kinetics

Michaelis–Menten and inhibition models integrate via `env.kinetics()` when bio rate laws are present:

```python
from ChemCompute import single_substrate_mm, competitive_inhibition

env = single_substrate_mm(s0=1.0, Vmax=1e-5, Km=1e-4)
checkpoints = env.kinetics(time=100.0, accuracy=0.01)

env = competitive_inhibition(s0=1.0, i0=0.2, Ki=1e-4)
```

Templates: `single_substrate_mm`, `competitive_inhibition`, `uncompetitive_inhibition`, `noncompetitive_inhibition`, `mixed_inhibition`, `sequential_pathway`.

---

## Testing

```bash
pytest tests/
python tests/manual/manual_validation.py   # 20 named equilibrium/kinetics cases
```

See also `tests/test_environment.py` for equilibrium, composition, buffering, titration, half-reactions, Pourbaix, activity, UV–Vis, and bio kinetics.

Kinetic plots from manual validation are written to `manual_test_output/kinetics/`.

---

## Half-reactions and Pourbaix diagrams

Half-reactions use **`@e`** as the electron token (never `e` or `e-` in `Reaction` strings). Species and concentrations mirror `Reaction`:

```python
from ChemCompute import Compound, Enviroment, HalfReaction, Pourbaix

hr = HalfReaction.from_string_simple_syntax(
    "Fe+3 + @e = Fe+2",
    concentrations=[0.01, 0.001],  # [Fe+3], [Fe+2] — no slot for @e
    E0=0.771,
)

env = Enviroment(
    hr,
    concentrations={"H+": 1e-7},
    buffer=["H+"],
)

env.set_electrode_potential(Eh=0.44)  # fixed E (V vs SHE)
env.equilibrium()

E = hr.E_at(env)  # Nernst E from concentrations

diagram = Pourbaix(env, pH_steps=30, Eh_steps=30).run()
diagram.plot(save="pourbaix.png", show=False)
print(diagram.junction_points[:3])  # species-labeled coordinates in model mode
```

### Speciation methods

| Method | Speed | Saved geometry | Use when |
|--------|-------|----------------|----------|
| `model` (default) | Fast | Analytic boundaries + junction points (`geometry_source='analytic'`) | Connected redox ladder per element, known E°/K/pKa, fixed `element_totals` |
| `equilibrium` | Slow | Dominance grid matrix (`geometry_source='grid'`) | Full coupling, stiff networks, or reactions not parsed by the graph model |

**Grid steps:** `pH_steps` / `Eh_steps` control the colored region grid in **both** modes (finer = less blocky fill). Boundary lines in `model` mode use `geometry_pH_steps` separately. **`progress=True`** prints scan percentage.

**Junction points:** numbered `P1`, `P2`, … — dominant-region coords via `diagram.junction_table()` (default) or full analytic set via `junction_table(source="analytic")`. Plot default: `boundary_mode="dominant"` (lines between neighboring regions only); use `boundary_mode="all"` for every analytic boundary.

**Model limits:** ideal dilute Nernst + pKa; independent per-element chains; oligomer/Ksp regions depend on totals and parsed reaction patterns; no cross-element redox.

```python
from ChemCompute import Pourbaix, build_pourbaix_graph

graph = build_pourbaix_graph(env)
diagram = Pourbaix(env, element_totals={"Se": 1.0}, speciation_method="model").run()
for junction in diagram.junction_points:
    print(junction.label, junction.pH, junction.Eh)
```

Complex syntax uses `=` and `&`:

```python
HalfReaction.from_string_complex_syntax(
    "Fe(OH)3.s & 3_H+ + @e = Fe+2 & 3_H2O.l",
    concentrations=[1.0, 1e-7, 0.05, 1.0],
    E0=-0.55,
)
```

Pass half-reactions directly to `Enviroment(rxn1, hr1, hr2, ...)`. With **two or more** half-reactions and no imposed `electrode_Eh`, equilibrium solves a **shared electrode potential** jointly with concentrations.

**Phase rules:** only explicit `.s` / `.l` solids and liquids are omitted from Q (activity 1). Undetermined phase stays in Q at concentration.

---

## Project layout

```
src/ChemCompute/
  _general.py      Compound, Reaction, Enviroment
  _equilibrium.py  Equilibrium solver and EquilibriumResult
  _kinetics.py     Time integration and plotting
  _activity.py     Ionic activity coefficients
  _buffer.py       Buffer capacity and Henderson-Hasselbalch
  _buffering.py    Solver-side constant-pH / species buffering
  _mixing.py       Environment combine and ScaledEnviroment
  _titration.py    Titration curves (sample + titrant environments)
  _half_reaction.py HalfReaction, Nernst electrode potential coupling
  _pourbaix.py     Pourbaix diagram scanner and PourbaixResult
  _uvvis.py        Beer-Lambert spectra
  _bio_kinetics.py Michaelis-Menten and inhibition integrator
  bio_templates/   Premade enzyme-kinetics environments
tests/
  helpers.py       Shared environment builders for tests
  test_general.py  Compound and Reaction
  test_environment.py  Enviroment API and features
  manual/          Reference environments and validation scripts
docs/index.md      Extended documentation
```

---

## License

MIT — see [LICENSE](LICENSE).

## Author

Mohammad Keifari — [mohammadkeifari2007@gmail.com](mailto:mohammadkeifari2007@gmail.com)
