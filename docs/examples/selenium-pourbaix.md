# Selenium Pourbaix diagram

Selenium at 1 M total Se, pH 0–10, with the aqueous acid–base ladder and the HSeO₄⁻ / H₂SeO₃ / Se / H₂Se redox chain.

## Code

```python
from ChemCompute import Enviroment, HalfReaction, Pourbaix, Reaction, XS

def ka(acid, base, pka):
    return Reaction.from_string(f"{acid} > H+ & {base}", K=10 ** (-pka))

env = Enviroment(
    Reaction.from_string("H2O.l > H+ & OH-", concentrations={"H2O": XS(0.0)}, K=1e-14),
    ka("HSeO4-", "SeO4-2", 1.92),
    ka("H2SeO3", "HSeO3-", 2.62),
    ka("HSeO3-", "SeO3-2", 7.19),
    ka("H2Se", "HSe-", 3.89),
    HalfReaction.from_string(
        "HSeO4- & 3_H+ & 2_@e = H2SeO3 & H2O.l", E0=1.15, name="HSeO4-/H2SeO3"
    ),
    HalfReaction.from_string(
        "H2SeO3 & 4_H+ & 4_@e = Se.s & 3_H2O.l", E0=0.74, name="H2SeO3/Se"
    ),
    HalfReaction.from_string("Se.s & 2_H+ & 2_@e = H2Se", E0=-0.11, name="Se/H2Se"),
    concentrations={"HSeO4-": 1.0},
    buffer=["H+"],
)

diagram = Pourbaix(
    env,
    pH_min=0.0,
    pH_max=10.0,
    Eh_min=-1.0,
    Eh_max=1.4,
    pH_steps=500,
    Eh_steps=500,
).run()

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

## Plot styles (in order)

1. **Filled** — colors each predominance region (default).
2. **Frame intersections** — filled regions plus markers where dominant-region boundaries meet the pH/Eh window.
3. **Predominance** — fill without junction markers.
4. **Boundaries** — line-only view with a compact legend.
5. **Labeled** — boundaries on white with dominant species written in each region.
6. **All boundaries** (`boundary_mode="all"`) — every analytic line, not only borders between neighboring regions.

## Figures

![Selenium Pourbaix, filled](../images/selenium_filled.png)

![Selenium Pourbaix with frame intersections](../images/selenium_frame.png)

![Selenium Pourbaix, predominance](../images/selenium_predominance.png)

![Selenium Pourbaix, boundaries](../images/selenium_boundaries.png)

![Selenium Pourbaix, labeled](../images/selenium_labeled.png)

![Selenium Pourbaix, all analytic boundaries](../images/selenium_all_boundaries.png)

## Related

- [Half-reactions and Pourbaix](../guide/pourbaix.md)
- [Half-reactions library](../libraries/half-reactions.md)
