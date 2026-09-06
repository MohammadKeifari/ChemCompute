# Half-reactions and Pourbaix diagrams

Half-reactions use **`@e`** as the electron token (never `e` or `e-` in `Reaction` strings).

## Half-reaction basics

```python
from ChemCompute import Compound, Enviroment, HalfReaction, Pourbaix

hr = HalfReaction.from_string(
    "Fe+3 & @e = Fe+2",
    concentrations=[0.01, 0.001],
    E0=0.771,
)

env = Enviroment(
    hr,
    concentrations={"H+": 1e-7},
    buffer=["H+"],
)

env.set_electrode_potential(Eh=0.44)  # fixed E (V vs SHE)
env.equilibrium()

E = hr.E_at(env)  # Nernst E from concentrations
```

Pass half-reactions directly to `Enviroment(rxn1, hr1, hr2, ...)`. With **two or more** half-reactions and no imposed `electrode_Eh`, equilibrium solves a **shared electrode potential** jointly with concentrations.

## Pourbaix scanner

```python
diagram = Pourbaix(env, pH_steps=30, Eh_steps=30).run()
diagram.plot(save="pourbaix.png", show=False)
print(diagram.junction_points[:3])
```

### Plot styles (recommended order)

```python
diagram.plot(plot_style="filled", save="selenium_filled.png", show=False)
diagram.plot(
    plot_style="filled",
    show_frame_intersections=True,
    save="selenium_frame.png",
    show=False,
)
diagram.plot_predominance(save="selenium_predominance.png", show=False)
diagram.plot_boundaries(save="selenium_boundaries.png", show=False)
diagram.plot(plot_style="labeled", save="selenium_labeled.png", show=False)
diagram.plot(
    plot_style="labeled",
    boundary_mode="all",
    save="selenium_all_boundaries.png",
    show=False,
)
```

| Style | Description |
|-------|-------------|
| `filled` | Colored predominance regions (default) |
| `filled` + `show_frame_intersections` | Filled + frame edge markers |
| `plot_predominance` | Fill without junction markers |
| `plot_boundaries` | Line-only with couple legend |
| `labeled` | White background + dominant species labels |
| `labeled` + `boundary_mode="all"` | Every analytic boundary |

For publication-quality fill, use `pH_steps=500`, `Eh_steps=500` (slower).

## Speciation methods

| Method | Speed | Saved geometry | Use when |
|--------|-------|----------------|----------|
| `model` (default) | Fast | Analytic boundaries + junction points | Connected redox ladder, known E°/K/pKa, fixed totals |
| `equilibrium` | Slow | Dominance grid matrix | Full coupling, stiff networks, unparsed reactions |

**Grid steps:**

| Parameter | `model` | `equilibrium` |
|-----------|---------|---------------|
| `pH_steps`, `Eh_steps` | Region fill/labels; finer = smoother | Layout **and** boundaries; finer = slower |
| `geometry_pH_steps` | Analytic line sampling (default 200) | N/A |

## Junction points

Numbered `P1`, `P2`, … — coordinates via `diagram.junction_table()` (dominant regions) or `junction_table(source="analytic")`.

Plot default: `boundary_mode="dominant"` (smooth clipped lines between neighboring regions). Use `boundary_mode="all"` for every analytic boundary.

## Graph model

```python
from ChemCompute import Pourbaix, build_pourbaix_graph

graph = build_pourbaix_graph(env)
diagram = Pourbaix(env, element_totals={"Se": 1.0}, speciation_method="model").run()
for junction in diagram.junction_points:
    print(junction.label, junction.pH, junction.Eh)
```

## Model limits

Ideal dilute Nernst + pKa; independent per-element chains; oligomer/Ksp regions depend on totals and parsed patterns; no cross-element redox in the graph model.

## Complex half-reaction syntax

```python
HalfReaction.from_string(
    "Fe(OH)3.s & 3_H+ & @e = Fe+2 & 3_H2O.l",
    concentrations=[1.0, 1e-7, 0.05, 1.0],
    E0=-0.55,
)
```

**Phase rules:** only explicit `.s` / `.l` solids and liquids are omitted from Q (activity 1).

## Library half-reactions

```python
from ChemCompute.half_reactions import hydrogen, oxygen, iron_iii, permanganate
```

See [Half-reactions library](../libraries/half-reactions.md) and [Selenium example](../examples/selenium-pourbaix.md).

## Related

- [Reaction syntax](syntax.md)
- [Equilibrium](../core/equilibrium.md)
