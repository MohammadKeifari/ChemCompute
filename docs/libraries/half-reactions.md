# Half-reactions library

`ChemCompute.half_reactions` names are functions that return shared `HalfReaction` objects with tabulated E° vs SHE (25 °C) and **all concentrations 0**.

## Usage

```python
from ChemCompute import Enviroment
from ChemCompute.half_reactions import hydrogen, oxygen, iron_iii, permanganate

env = Enviroment(
    hydrogen(),
    oxygen(),
    iron_iii(),
    concentrations={"Fe+3": 0.01, "Fe+2": 0.001},
    buffer=["H+"],
)
```

## Categories

| Module | Examples |
|--------|----------|
| `_water.py` | `hydrogen` (SHE, E°=0), `she`, `oxygen`, `oxygen_hydroxide`, peroxide couples |
| `_halogens.py` | `fluorine`, `chlorine`, `bromine`, `iodine`, `triiodide`, `hypochlorite` |
| `_metals.py` | `iron_iii`, `iron`, `copper`, `silver`, `zinc`, `aluminum`, … |
| `_oxoanions.py` | `permanganate`, `dichromate`, `nitrate`, `sulfate`, … |
| `_complexes.py` | `ferricyanide` — Fe(CN)6³⁻/Fe(CN)6⁴⁻ |

## Conventions

- Library ions and liquid water via `.token`
- Metal solids written as `M.s` (not library compounds)
- Aqueous peroxide/sulfur dioxide as `H2O2.aq`, `SO2.aq` so they stay in Q
- `hydrogen()` is the standard hydrogen electrode; `she` is the same function

## API

```python
from ChemCompute.half_reactions import all_half_reactions
```

## Related

- [Half-reactions and Pourbaix](../guide/pourbaix.md)
- [Selenium example](../examples/selenium-pourbaix.md)
