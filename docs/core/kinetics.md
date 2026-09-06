# Kinetics

`env.kinetics()` integrates mass-action rate laws forward in time from the current concentrations and rate constants `kf` / `kb`.

## Basic usage

```python
results = env.kinetics(
    time=10.0,
    accuracy=1e-3,
    checkpoint_time=[1.0, 5.0, 10.0],
)
final = results[-1]
```

`accuracy` is the integration time step.

## Plotting

```python
env.kinetics(
    time=10.0,
    plot="interactive",  # or "save" or False
    directory="./plot.png",
    colors=["#26547c", "#ef476f"],
)
```

One color per compound when `colors` is provided.

## From equilibrium

Kinetics uses current `env.concentrations`. Run equilibrium first to integrate from an equilibrated state:

```python
env.apply_equilibrium(method="newton", tol=1e-10)
env.kinetics(time=5.0, plot="save", directory="approach.png")
```

## Buffered species

Species in `env.buffer` stay fixed during integration (same mechanism as equilibrium buffering).

## Bio kinetics

When bio rate laws are present, `env.kinetics()` delegates to the bio integrator. See [Bio kinetics](../guide/bio-kinetics.md).

## Example: Jungle Model

Four irreversible steps with autocatalytic growth:

```python
env = Enviroment(
    Reaction.from_string("R & T > 2_R", kf=0.01, kb=0.0, K=1e12),
    Reaction.from_string("R > D", kf=0.5, kb=0.0, K=1e12),
    Reaction.from_string("R & W > 2_W", kf=0.01, kb=0.0, K=1e12),
    Reaction.from_string("W > inert", kf=0.5, kb=0.0, K=1e12),
    concentrations={"R": 70.0, "T": 100.0, "W": 20.0},
)
env.kinetics(time=40.0, accuracy=0.05, plot="save", directory="jungle_model.png")
```

Full walkthrough: [Jungle Model example](../examples/jungle-model.md).

## Related

- [Reaction](reaction.md) — rate constants and temperature dependence
- [Temperature-dependent K and k](../guide/syntax.md)
