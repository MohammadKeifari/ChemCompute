# Compounds library

`ChemCompute.compounds` provides common molecules, ions, coordination complexes, and sparingly soluble salts.

## Usage

```python
from ChemCompute.compounds import water, h_plus, oh_minus, mno4, fescn, agcl, get

rxn = Reaction.from_string(
    f"{water().token} > {h_plus().token} & {oh_minus().token}",
    concentrations=[1.0, 1e-7, 1e-7],
    K=1e-14,
)
assert rxn.reactants[0]["compound"] is water()
assert get("H2O") is water()
```

## Categories

| Module | Contents |
|--------|----------|
| `_molecules.py` | Diatomics, solvents, acids, organics — `h2o`/`water`, `ch3cooh`, `nh3`, … |
| `_ions.py` | Aqueous ions and colored oxoanions — `h_plus`, `mno4`, `cro4`, `i3`, … |
| `_complexes.py` | FeSCN²⁺, ammines, ferroin, hexacyanoferrates, dyes — `fescn`, `ferroin`, … |
| `_solids.py` | Ksp salts — `agcl`, `caf2`, `feoh3`, … |

## Physical properties

- `mp` and `bp` in **kelvin** (1 atm)
- Ions and decomposing species often omit mp/bp

## UV-Vis

Connected envelope spectra where a clear aqueous trace exists (`mno4`, `fescn`, `ferroin`, …). ε in M⁻¹ m⁻¹; default 0.01 m path = 1 cm cuvette.

See [UV-Vis guide](../guide/uvvis.md).

## Aliases

Some names alias the same object (`water` / `h2o`, `ethanol` / `c2h5oh`).

## API

```python
all_compounds()  # unique Compound instances
get("SO4-2")     # lookup by formula
```

## Related

- [Compound](../core/compound.md)
- [Reactions library](reactions.md)
