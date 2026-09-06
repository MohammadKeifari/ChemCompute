# Installation

## From PyPI

```bash
pip install chemcompute
```

This installs NumPy and Matplotlib as dependencies.

## From source

```bash
git clone https://github.com/MohammadKeifari/ChemCompute.git
cd ChemCompute
pip install -e .
```

## Development and documentation

For tests and local documentation preview:

```bash
pip install -e ".[dev,docs]"
pytest tests/
mkdocs serve
```

Open `http://127.0.0.1:8000` for the documentation site.

## Requirements

- Python **3.9+**
- NumPy ≥ 1.19
- Matplotlib ≥ 3.3 (plotting: kinetics, titration, Pourbaix, UV-Vis)

## Import

```python
from ChemCompute import Compound, Reaction, Enviroment, EquilibriumResult
from ChemCompute import HalfReaction, Pourbaix, Titration, XS
```

`Environment` is an alias for `Enviroment`.

## Next steps

Continue with the [Quick start](quickstart.md) tutorial.
