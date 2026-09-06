# Ag₂T precipitation titration

50 mL of 0.07 M AgF titrated with 0.01 M (NH₄)₂T. Precipitation of Ag₂T (Ksp = 4×10⁻¹²) is predicted when Qsp = [Ag⁺]² [T²⁻] crosses Ksp, near **0.0042 µL** titrant added.

## Code

```python
from ChemCompute import Enviroment, Reaction, Titration, XS
import numpy as np
import matplotlib.pyplot as plt

agf = Enviroment(
    Reaction.from_string(
        "H2O.l > H+ & OH-",
        concentrations={"H2O": XS(55.5), "H+": 1e-7, "OH-": 1e-7},
        K=1e-14,
    ),
    Reaction.from_string("HF.aq > H+ & F-", K=10 ** (-3.1)),
    concentrations={"Ag+": 0.07, "F-": 0.07},
    volume=0.050,
)
titrant = Enviroment(
    Reaction.from_string("H2T > H+ & HT-", K=10 ** (-4.2)),
    Reaction.from_string("HT- > H+ & T-2", K=10 ** (-6.5)),
    Reaction.from_string("NH4+ > H+ & NH3", K=10 ** (-9.2)),
    Reaction.from_string("Ag+ & 2_NH3 > Ag(NH3)2+", K=2e7),
    concentrations={"NH4+": 0.02, "T-2": 0.01},
    volume=1.0,
)

volumes_ul = np.linspace(0.0, 0.012, 61)
curve = Titration(agf, titrant, volumes=volumes_ul * 1e-6).run(
    method="newton", tol=1e-10
)
ksp = 4e-12
qsp = curve.species("Ag+") ** 2 * curve.species("T-2")

# Custom Qsp plot (micro-liter scale)
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(volumes_ul, qsp, label=r"$Q_\mathrm{sp}=[\mathrm{Ag}^+]^2[\mathrm{T}^{2-}]$")
ax.axhline(ksp, linestyle="--", label=r"$K_\mathrm{sp}=4\times10^{-12}$")
ax.set_xlabel("Titrant volume added (µL)")
ax.set_ylabel(r"$Q_\mathrm{sp}$")
ax.set_yscale("log")
ax.legend()
fig.savefig("ag2t_titration_qsp.png", dpi=150)
```

The figure below uses the same data as `docs/generate_readme_figures.py`.

## Figure

![Ag₂T precipitation onset from Qsp vs Ksp](../images/ag2t_titration_qsp.png)

## Related

- [Mixing and titration](../guide/mixing-and-titration.md)
- [Equilibrium](../core/equilibrium.md)
