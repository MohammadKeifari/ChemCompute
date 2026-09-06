# Troubleshooting

## Equilibrium solver fails or stalls

**Symptoms:** `RuntimeError`, non-convergence, or concentrations stuck at `min_concentration`.

**Checks:**

1. **Stoichiometry** — reactant and product formulas must balance. Use `_` for stoichiometric coefficients (`2_H2O`, not `2H2O` in the token grammar).
2. **Initial guess** — set reasonable starting concentrations on the `Enviroment` or individual reactions.
3. **Method and tolerance** — default is `method="newton"` with `tol=1e-10`. Try lowering tolerance slightly or increasing `max_iter`.
4. **Strong acids/bases** — use `infinite_K=True` for full dissociation instead of a finite K that fights the solver.
5. **XS for excess solids/liquids** — mark bulk water or excess precipitate with `XS(...)` so activity stays 1 without consuming the variable.

```python
result = env.equilibrium(method="newton", tol=1e-10, max_iter=8000, return_details=True)
print(result.message)
```

## Stoichiometry or parse errors

**Symptoms:** `ValueError` when calling `Reaction.from_string` or `HalfReaction.from_string`.

**Checks:**

- Use `>` for reactions (not `=` except in half-reactions).
- Separate species with `&` on each side.
- Phase suffixes: `.aq`, `.s`, `.l`, `.g` — see [Reaction syntax](guide/syntax.md).
- Half-reactions: `@e` on the **left** only; never `e-` in reaction strings.

## Kinetics: species stay at zero

**Symptoms:** A species never appears in kinetic plots despite being in the mechanism.

**Cause:** Initial concentrations were set only on a reaction string, not on the environment.

**Fix:**

```python
env = Enviroment(
    Reaction.from_string("R & T > 2_R", kf=0.01, kb=0.0, K=1e12),
    concentrations={"R": 70.0, "T": 100.0, "W": 20.0},  # set here
)
```

## Pourbaix diagram looks wrong

1. **Buffer H⁺** — pass `buffer=["H+"]` when pH should float with Nernst lines.
2. **`model` vs `equilibrium`** — default `model` is fast but assumes parsed redox ladders; use `speciation_method="equilibrium"` for fully coupled networks (slower).
3. **Grid resolution** — increase `pH_steps` and `Eh_steps` for smoother fills; analytic boundaries use `geometry_pH_steps` separately in `model` mode.
4. **Element totals** — pass `element_totals={"Se": 1.0}` when totals are not inferred from concentrations.

See [Pourbaix guide](guide/pourbaix.md).

## Titration volume scale

For very small titrant volumes (µL range), the default `Titration.plot` x-axis may be hard to read. Extract species arrays and plot Qsp or pH vs volume manually — see [Ag₂T example](examples/ag2t-titration.md).

## Library object confusion

Library names are **functions**: call `water()`, not `water`. Concentrations belong on the environment, not baked into library singletons.

## Still stuck?

Open an issue on [GitHub](https://github.com/MohammadKeifari/ChemCompute/issues) with a minimal reproducible `Enviroment` definition.
