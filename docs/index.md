# ChemCompute documentation

ChemCompute models multi-reaction chemical systems in Python. You define compounds and reactions once, wrap them in an `Enviroment`, and then run:

- **Equilibrium** — solve for concentrations where each reaction satisfies its mass-action expression (Q/K)
- **Kinetics** — integrate concentrations forward in time from rate laws
- **Titration** — mix a sample with a titrant over a volume grid and equilibrate at each step
- **Pourbaix** — map predominance vs pH and electrode potential

Both equilibrium and kinetics share the same reaction network, stoichiometry, and concentration state.

## Where to start

| If you want to… | Start here |
|-----------------|------------|
| Install the package | [Installation](getting-started/installation.md) |
| Run your first calculation | [Quick start](getting-started/quickstart.md) |
| Understand the object model | [Concepts](concepts.md) |
| See full worked examples with figures | [Examples](examples/index.md) or the [README gallery](https://github.com/MohammadKeifari/ChemCompute#examples) |

## Core components

- [Compound](core/compound.md) — formulas, phases, melting/boiling points, UV-Vis spectra
- [Reaction](core/reaction.md) — stoichiometry, K, rate constants, string notation
- [Environment](core/environment.md) — multi-reaction systems, mixing, buffering
- [Equilibrium](core/equilibrium.md) — Newton / BGD / SGD solvers, phases, excess
- [Kinetics](core/kinetics.md) — time integration and plotting

## Guides

- [Reaction syntax](guide/syntax.md) — `&`, `@e`, phase suffixes, live compound tokens
- [Mixing and titration](guide/mixing-and-titration.md)
- [Activity and buffering](guide/activity-and-buffering.md)
- [Half-reactions and Pourbaix](guide/pourbaix.md)
- [UV-Vis](guide/uvvis.md)
- [Bio kinetics](guide/bio-kinetics.md)

## Libraries

Pre-built compounds, reactions, half-reactions, and coupled environments with tabulated K or E° and zero concentrations:

- [Libraries overview](libraries/overview.md)
- [Compounds](libraries/compounds.md)
- [Reactions](libraries/reactions.md)
- [Half-reactions](libraries/half-reactions.md)
- [Environments](libraries/environments.md)

## Requirements

- Python 3.9 or higher
- NumPy
- Matplotlib (for plotting)

## License

MIT — see [LICENSE](https://github.com/MohammadKeifari/ChemCompute/blob/main/LICENSE).
