# Bio kinetics

Michaelis–Menten and inhibition models integrate through `env.kinetics()` when bio rate laws are present on reactions.

## Templates

```python
from ChemCompute import single_substrate_mm, competitive_inhibition

env = single_substrate_mm(s0=1.0, Vmax=1e-5, Km=1e-4)
checkpoints = env.kinetics(time=100.0, accuracy=0.01)

env = competitive_inhibition(s0=1.0, i0=0.2, Ki=1e-4)
```

Available templates in `ChemCompute.bio_templates`:

- `single_substrate_mm`
- `competitive_inhibition`
- `uncompetitive_inhibition`
- `noncompetitive_inhibition`
- `mixed_inhibition`
- `sequential_pathway`

## Low-level API

```python
from ChemCompute import integrate_bio_kinetics, uses_bio_kinetics, RATE_LAW_FUNCTIONS
```

`uses_bio_kinetics(env)` returns whether the environment will use the bio integrator.

## Related

- [Kinetics](../core/kinetics.md)
- [Project layout](../project/layout.md) — `bio_templates/` module
