# Project layout

```
src/ChemCompute/
  _general.py       Compound, Reaction, Enviroment
  _equilibrium.py   Equilibrium solver and EquilibriumResult
  _kinetics.py      Time integration and plotting
  _activity.py      Ionic activity coefficients
  _buffer.py        Buffer capacity and Henderson-Hasselbalch
  _buffering.py     Solver-side constant-pH / species buffering
  _mixing.py        Environment combine and ScaledEnviroment
  _titration.py     Titration curves (sample + titrant environments)
  _half_reaction.py HalfReaction, Nernst electrode potential coupling
  _pourbaix.py      Pourbaix diagram scanner and PourbaixResult
  _uvvis.py         Beer-Lambert spectra
  _bio_kinetics.py  Michaelis-Menten and inhibition integrator
  bio_templates/    Premade enzyme-kinetics environments
  compounds/        Library of common molecules, ions, and complexes
  reactions/        Library of acid–base, precipitation, complex, and redox reactions
  half_reactions/   Library of aqueous half-reactions vs SHE
  environments/     Coupled polyprotic, complex, precipitation, and redox environments
tests/
  helpers.py        Shared environment builders for tests
  test_general.py   Compound and Reaction
  test_environment.py  Enviroment API and features
  manual/           Reference environments and validation scripts
docs/               MkDocs source (this site)
  images/           README and docs gallery figures
  generate_readme_figures.py
mkdocs.yml          Site configuration
```

## Package entry point

Public API is re-exported from `ChemCompute/__init__.py`. Library subpackages (`compounds`, `reactions`, `half_reactions`, `environments`) expose cached getter functions.

## Related

- [Testing](testing.md)
- [Contributing](../contributing.md)
