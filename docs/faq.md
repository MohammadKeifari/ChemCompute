# FAQ

## Why is it spelled `Enviroment`?

Historical naming in the codebase. `Environment` is exported as an alias:

```python
from ChemCompute import Enviroment, Environment
assert Enviroment is Environment
```

## `@e` vs `e-` in strings

Electrons in **half-reactions** use the token `@e` on the left side of the equation:

```python
HalfReaction.from_string("Fe+3 & @e = Fe+2", E0=0.771)
```

Do **not** write `e-` in `Reaction` or `HalfReaction` strings — the parser treats `@e` specially and excludes it from concentration slots.

## When should I use Pourbaix `model` vs `equilibrium`?

| Method | When to use |
|--------|-------------|
| `model` (default) | Known E°, pKa, and redox ladder per element; fast analytic boundaries |
| `equilibrium` | Full coupling, stiff networks, or chemistry the graph model cannot parse |

See [Pourbaix guide](guide/pourbaix.md).

## What is the default equilibrium method?

`newton` with tolerance `1e-10`. Older docs may mention other methods; the current default is Newton-based root finding on mass-action residuals.

## How do library getters work?

Each name (`water`, `iron_iii`, `water_kw`, …) is a function returning a **cached singleton**. Call it every time:

```python
from ChemCompute.compounds import water
rxn = Reaction.from_string(f"{water().token} > H+ & OH-", K=1e-14)
```

Concentrations start at zero on library objects; set amounts on the `Enviroment`.

## What does `XS` mean?

**Excess** — marks a solid or liquid present in bulk so its activity is fixed at 1 without being a free concentration variable:

```python
Reaction.from_string("CaF2.s > Ca+2 & 2_F-", {"CaF2": XS(10.0)}, K=5e-9)
```

## Can I mix two solutions?

Yes — `combine_environments` or `Titration` for sequential addition. See [Mixing and titration](guide/mixing-and-titration.md).

## Where are the README gallery figures?

Under `docs/images/`. Regenerate with:

```bash
python docs/generate_readme_figures.py
```

## Related

- [Troubleshooting](troubleshooting.md)
- [Concepts](concepts.md)
