# Contributing

Thank you for improving ChemCompute. This page covers local setup, tests, documentation, and figure regeneration.

## Development install

```bash
git clone https://github.com/MohammadKeifari/ChemCompute.git
cd ChemCompute
pip install -e ".[dev,docs]"
```

## Run tests

```bash
pytest tests/
python tests/manual/manual_validation.py
```

See [Testing](project/testing.md) for module-level coverage.

## Documentation site

Documentation uses [MkDocs](https://www.mkdocs.org/) with the Material theme. Source files live in `docs/`; configuration is in `mkdocs.yml` at the repository root.

### Local preview

```bash
mkdocs serve
```

Open `http://127.0.0.1:8000`.

### Strict build

Before opening a pull request that changes docs:

```bash
mkdocs build --strict
```

This fails on broken internal links or nav issues. Build output goes to `site/` (gitignored).

### Adding pages

1. Create a Markdown file under `docs/`.
2. Add an entry to the `nav:` section in `mkdocs.yml`.
3. Run `mkdocs build --strict`.

Keep the [README](https://github.com/MohammadKeifari/ChemCompute/blob/main/README.md) as a visual showcase; put detailed prose in `docs/`.

## Regenerate gallery figures

The README and [Examples](examples/index.md) embed PNGs from `docs/images/`:

```bash
python docs/generate_readme_figures.py
```

Requires Matplotlib and a working ChemCompute install (the script adds `tests/manual` to `sys.path` for `build_env16`).

## Code conventions

- Match existing style in surrounding modules.
- Library compounds/reactions/half-reactions: cached getters, zero concentrations, `.token` in strings.
- Half-reactions: `@e` token, E° vs SHE at 25 °C unless noted.
- Run `pytest tests/` for changes touching solvers or parsing.

## License

Contributions are under the project MIT license — see [LICENSE](https://github.com/MohammadKeifari/ChemCompute/blob/main/LICENSE).

## Related

- [Project layout](project/layout.md)
- [FAQ](faq.md)
