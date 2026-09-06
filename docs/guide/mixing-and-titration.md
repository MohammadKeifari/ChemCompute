# Mixing and titration

## Mixing environments

Volume-weighted mixing combines concentrations and reaction networks:

```python
envA = Enviroment(rxn_a, concentrations={"A": 1.0}, volume=1.0)
envB = Enviroment.from_compounds({"B": 0.1}, volume=1.0)

envC = envA + envB                         # volume = 2.0
envD = 0.5 * envA + 4 * envB               # effective volume = 4.5 L
envE = Enviroment.combine((0.5, envA), (4, envB))

envF = envC.add_compounds({"A": 1.0}, volume=1.0)
envG = envC.add_compounds({"A": 1.0}, volume=1.0, coefficient=4.0)
```

Mixing rule: `effective_volume = coeff × volume`, total moles per species are conserved, final concentration = moles / total effective volume.

`ScaledEnviroment` wraps `(coefficient, env)` for explicit combine calls.

## Titration

Mix a **sample** with a **titrant** over titrant volumes, equilibrate at each step:

```python
from ChemCompute import Titration, Enviroment

sample = Enviroment(...)  # set sample.volume (e.g. 0.050 L)
titrant = Enviroment(...)  # stock concentrations; volume defines reference size

curve = Titration(
    sample,
    titrant,
    volume_min=0.0,
    volume_max=0.05,
    steps=100,
).run(method="newton", tol=1e-10)
```

Or pass explicit volumes (e.g. microlitres as litres):

```python
import numpy as np
volumes_ul = np.linspace(0.0, 0.012, 61)
curve = Titration(sample, titrant, volumes=volumes_ul * 1e-6).run(method="newton", tol=1e-10)
```

### Accessing results

```python
curve.titrant_volumes     # L added at each step
curve.pH
curve.matrix()            # shape (n_steps, n_compounds)
curve.species("Ag+")      # one species vs volume
curve.speciation          # list of {formula: c} per step
```

### Plotting

```python
curve.plot(species=["H+", "OH-"], plot="save", directory="titration.png")
curve.plot_pH(plot="interactive")
```

The sample and titrant are not mutated.

## Example

[Ag2T precipitation titration](../examples/ag2t-titration.md) — 50 mL 0.07 M AgF with 0.01 M (NH4)2T.

## Related

- [Environment](../core/environment.md)
- [Equilibrium](../core/equilibrium.md)
