# Environments library

`ChemCompute.environments` couples related library reactions (and sometimes half-reactions) into a **new** `Enviroment` each call. Concentrations start at 0 unless you pass `concentrations=`.

## Usage

```python
from ChemCompute.environments import (
    phosphoric_acid,
    silver_chloride,
    copper_hydroxide_ammine,
    water_limits,
)

h3po4 = phosphoric_acid(concentrations={"H3PO4": 0.10})
agcl = silver_chloride(concentrations={"Ag+": 1e-3, "Cl-": 0.10})
window = water_limits()
```

Each environment uses the same shared library objects (`water()`, `water_kw()`, `agcl()`, …). Amounts live on that environment instance.

## Acid–base (`_acids.py`)

Polyprotic acids with every deprotonation plus Kw:

- `phosphoric_acid`, `carbonic_acid`, `hydrogen_sulfide`, `sulfurous_acid`, `sulfuric_acid`, `ammonia`

## Precipitation and speciation (`_precipitation.py`)

- `silver_chloride` — AgCl(s) Ksp + AgCl2⁻ formation
- `silver_chloride_ammonia` — adds ammine chemistry
- `copper_hydroxide_ammine`, `zinc_hydroxide_ammine`
- `calcium_carbonate`, `calcium_fluoride`, `barium_sulfate`, …

## Complexes (`_complexes.py`)

- `iron_thiocyanate`, `ferroin`, ammine networks, cyanide complexes, triiodide

## Redox (`_redox.py`)

- `water_limits` — O₂/H₂O + H⁺/H₂ window with Kw
- `iron_couple` — Fe³⁺/Fe²⁺ with Kw
- `daniel_cell` — Cu²⁺/Cu with Zn²⁺/Zn

## Factory signature

```python
def phosphoric_acid(concentrations=None, *, T=298, volume=1.0):
    ...
```

## API

```python
from ChemCompute.environments import all_environments
```

## Related

- [Environment](../core/environment.md)
- [Reactions library](reactions.md)
- [Half-reactions library](half-reactions.md)
