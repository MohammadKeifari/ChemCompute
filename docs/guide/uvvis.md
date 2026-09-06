# UV-Vis (Beer–Lambert)

Molar absorptivity is stored on the `Compound` as a `SpectrumSpec`.

## User-defined spectrum

```python
from ChemCompute import Compound, SpectrumSpec, uvvis_spectrum

inh = Compound(
    "InH",
    phase_point_list=[{"phase": "aq", "temperature": 298}],
    spectrum=SpectrumSpec(
        points=[(400e-9, 12000), (450e-9, 25000), (500e-9, 8000)],
        extrapolate="flat",  # or "none" for zero off-tabulated wavelengths
    ),
)
```

Attach after parsing:

```python
env.set_spectrum("InH", SpectrumSpec(...))
A = uvvis_spectrum(env, wavelengths=[450e-9, 500e-9], path_length=0.01)
```

## Units

- Wavelength in **metres**
- ε in **M⁻¹ m⁻¹** (literature M⁻¹ cm⁻¹ × 100)
- Default `path_length=0.01` m → 1 cm cuvette

Library spectra use connected (wavelength, ε) points — piecewise-linear envelopes, not single spikes.

## Library compounds

```python
from ChemCompute import Reaction, uvvis_spectrum
from ChemCompute.compounds import water, h_plus, oh_minus, mno4, fescn

rxn = Reaction.from_string(
    f"{water().token} > {h_plus().token} & {oh_minus().token}",
    concentrations=[1.0, 1e-7, 1e-7],
    K=1e-14,
)
A = uvvis_spectrum(env, wavelengths=[500e-9, 525e-9, 550e-9], path_length=0.01)
```

`mno4()` traces the permanganate visible band (peak ~525 nm). `fescn()` is the FeSCN²⁺ LMCT envelope (~447 nm).

Look up species with `compounds.get("H2O")`.

## Related

- [Compound](../core/compound.md)
- [Compounds library](../libraries/compounds.md)
