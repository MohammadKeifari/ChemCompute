# Compound

The `Compound` class represents a chemical species with formula, phase behaviour, and optional physical and spectroscopic data.

## Construction

```python
from ChemCompute import Compound

water = Compound(
    "H2O",
    phase_point_list=[{"phase": "l", "temperature": 298}],
    mp=273.15,
    bp=373.15,
)

h_plus = Compound("H+", phase_point_list=[{"phase": "aq", "temperature": 298}], charge=1)
```

## Phase

Call `compound.phase(T)` to resolve phase at temperature `T` (kelvin). Phases: `s`, `l`, `g`, `aq`.

In reaction strings, append `.s`, `.l`, `.g`, or `.aq` to the species token. Solids and liquids are omitted from equilibrium Q (activity 1).

## Melting and boiling points

`mp` and `bp` are 1 atm values in **kelvin** (same scale as `env.T`). Ions and decomposing species often omit them.

## Charge

Set `charge` explicitly or infer from trailing `+` / `-` in formulas (`H+`, `SO4-2`, `Fe(CN)6-4`). Non-zero charges populate `env.charge_map` when the compound enters an environment.

## UV-Vis spectrum

Attach a `SpectrumSpec` for Beer-Lambert calculations:

```python
from ChemCompute import Compound, SpectrumSpec

dye = Compound(
    "D",
    phase_point_list=[{"phase": "aq", "temperature": 298}],
    spectrum=SpectrumSpec(
        points=[(500e-9, 1000.0)],
        extrapolate="flat",
    ),
)
```

After building an environment from strings, attach spectra by formula:

```python
env.set_spectrum("FeSCN+2", SpectrumSpec(...))
```

## Live compound token

Library and custom compounds can be interpolated into reaction strings via `.token` so the same object (with spectrum, phase, mp/bp) is reused:

```python
from ChemCompute.compounds import water

rxn = Reaction.from_string(
    f"{water().token} > H+ & OH-",
    concentrations={"H2O": XS(55.5)},
    K=1e-14,
)
assert rxn.reactants[0]["compound"] is water()
```

Do **not** use `f"{water()}"` — that prints the formula, not the live slot.

## Library compounds

```python
from ChemCompute.compounds import water, h_plus, mno4, get

assert water() is water()  # same object every call
assert get("H2O") is water()
```

See [Compounds library](../libraries/compounds.md).

## Related

- [Reaction syntax](../guide/syntax.md) — phase suffixes in strings
- [UV-Vis guide](../guide/uvvis.md)
