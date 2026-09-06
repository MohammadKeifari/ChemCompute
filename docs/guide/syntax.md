# Reaction syntax

Use `Reaction.from_string(...)` and `HalfReaction.from_string(...)` with a single grammar.

## Token reference

| Role | Token | Example |
|------|-------|---------|
| Species on one side | `&` | `HSeO4- & 3_H+ & 2_@e` |
| Reaction direction | `>` | `HA > H+ & A-` |
| Half-reaction sides | `=` | `Ox = Red` |
| Stoichiometry | prefix `n_` | `3_H+`, `2_@e` |
| Rate order (Reaction) | suffix `_n` | `A_2` (optional; default = stoichiometry) |
| Phase | suffix | `.aq`, `.s`, `.l`, `.g` |
| Electrons | `@e` only | never bare `e-` in Reaction |
| Live compound | `{water().token}` | keeps spectrum, phase, mp/bp |

## Charge inference

Ionic charge is inferred from trailing `+` / `-` in species names (`H+`, `SeO4-2`, `Fe(CN)6-4`). Explicit `Compound(..., charge=...)` overrides when needed.

## Concentrations

```python
from ChemCompute import XS

kw = Reaction.from_string(
    "H2O.l > H+ & OH-",
    K=1e-14,
    concentrations={"H2O": XS(55.5), "H+": 0, "OH-": 0},
)
```

Legacy list form: reactants then products in order.

## Half-reactions

Electrons appear only on the **left** (oxidized) side:

```python
from ChemCompute import HalfReaction

hr = HalfReaction.from_string(
    "Fe+3 & @e = Fe+2",
    concentrations={"Fe+3": 0.01, "Fe+2": 0.001},
    E0=0.771,
)
```

Omit `@e` from concentration lists — only count non-electron species.

## Live library compounds

```python
from ChemCompute.compounds import water

kw = Reaction.from_string(
    f"{water().token} > H+ & OH-",
    concentrations={"H2O": XS(55.5)},
    K=1e-14,
)
assert kw.reactants[0]["compound"] is water()
```

## Migration from old syntax

Replace `+` between species with `&`:

- `A + B > C` → `A & B > C`
- `Fe+3 + @e = Fe+2` → `Fe+3 & @e = Fe+2`

## Temperature-dependent parameters

```python
rxn = Reaction.from_string(
    "A > B",
    K=2.0, kf=0.5, kb=0.25,
    enthalpy=-50000,
    activation_energy_forward=50000,
    activation_energy_backward=100000,
    T=298,
)
rxn.T = 350
```

When `entropy` is set, van't Hoff uses ΔG = ΔH − TΔS.

## Related

- [Compound](../core/compound.md)
- [Half-reactions and Pourbaix](pourbaix.md)
