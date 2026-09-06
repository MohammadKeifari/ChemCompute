# Reactions library

`ChemCompute.reactions` names are functions that return shared `Reaction` objects with tabulated K (or `infinite_K`) and **all concentrations 0**.

## Usage

```python
from ChemCompute import Enviroment
from ChemCompute.reactions import water_kw, acetic_acid, agcl_ksp, fescn_kf

env = Enviroment(
    water_kw(),
    acetic_acid(),
    agcl_ksp(),
    fescn_kf(),
    concentrations={"CH3COOH": 0.10, "Ag+": 1e-3},
)
```

## Categories

### Acid–base (`_acid_base.py`)

- `water_kw` — Kw = 1e-14
- Strong acids: `hydrochloric_acid`, `nitric_acid`, `perchloric_acid`, …
- Weak acids/bases: `acetic_acid`, `ammonia`, polyprotic steps (`phosphoric_acid_1`, …)

### Precipitation (`_precipitation.py`)

Ksp dissolutions for common salts: `agcl_ksp`, `caf2_ksp`, `feoh3_ksp`, `pbso4_ksp`, …

### Complexes (`_complexes.py`)

Overall formation: `fescn_kf`, `agcl2_kf`, ammines, `ferroin_kf`, triiodide, hexacyanoferrates, …

### Redox titrations (`_redox.py`)

Analytical `infinite_K` titrations: `permanganate_iron`, `dichromate_iron`, `iodine_thiosulfate`, …

## Phase conventions

- Liquid water from `water()` — activity 1, omitted from Q
- Solids from library compounds — activity 1
- Molecular acids in Q written as `.aq` when stored as gas/liquid in compounds

## API

```python
from ChemCompute.reactions import all_reactions

for rxn in all_reactions():
    ...
```

## Related

- [Reaction](../core/reaction.md)
- [Environments library](environments.md)
