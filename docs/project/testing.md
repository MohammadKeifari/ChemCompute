# Testing

## Automated tests

```bash
pytest tests/
```

Pytest configuration lives in `pyproject.toml` under `[tool.pytest.ini_options]`. Tests are under `tests/` with the `test_*.py` naming convention.

### Main test modules

| File | Coverage |
|------|----------|
| `test_general.py` | `Compound`, `Reaction`, string parsing |
| `test_environment.py` | Equilibrium, composition, buffering, titration, half-reactions, Pourbaix, activity, UV-Vis, bio kinetics |
| `test_half_reactions.py` | Half-reaction library getters |
| `test_pourbaix.py` | Pourbaix diagram behavior |

## Manual validation

Twenty named equilibrium and kinetics cases:

```bash
python tests/manual/manual_validation.py
```

Reference environment builders (`build_env16` … `build_env20`) are in `tests/manual/complex_environments.py`. Expected concentrations are checked against stored values in the manual suite.

Kinetic plots from manual validation are written to `manual_test_output/kinetics/`.

## Coverage (optional)

With dev dependencies installed:

```bash
pip install -e ".[dev]"
pytest tests/ --cov=ChemCompute --cov-report=term-missing
```

## Related

- [Env 16 example](../examples/env16-equilibrium.md)
- [Contributing](../contributing.md)
