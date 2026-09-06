# Jungle Model kinetics

Autocatalytic growth of R on resource T, death of R to D, conversion of R into W, and death of W to inert:

```
R + T  ->  2R       k = 0.01
R      ->  D        k = 0.5
R + W  ->  2W       k = 0.01
W      ->  inert    k = 0.5
```

## Code

```python
from ChemCompute import Enviroment, Reaction

env = Enviroment(
    Reaction.from_string("R & T > 2_R", kf=0.01, kb=0.0, K=1e12),
    Reaction.from_string("R > D", kf=0.5, kb=0.0, K=1e12),
    Reaction.from_string("R & W > 2_W", kf=0.01, kb=0.0, K=1e12),
    Reaction.from_string("W > inert", kf=0.5, kb=0.0, K=1e12),
    concentrations={"R": 70.0, "T": 100.0, "W": 20.0},
)
env.kinetics(
    time=40.0,
    accuracy=0.05,
    plot="save",
    directory="jungle_model.png",
    colors=["#2a9d8f", "#e9c46a", "#6d6875", "#e76f51", "#264653"],
)
```

Pass initial concentrations on the `Enviroment`, not only on individual reaction strings, so every species starts at the intended value.

## Figure

![Jungle Model kinetic trajectories](../images/jungle_model.png)

## Related

- [Kinetics](../core/kinetics.md)
- [Reaction syntax](../guide/syntax.md)
