"""Regenerate README gallery figures under docs/images/."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests" / "manual"))

from ChemCompute import (  # noqa: E402
    Compound,
    Enviroment,
    HalfReaction,
    Pourbaix,
    Reaction,
    Titration,
    XS,
)
from complex_environments import build_env16  # noqa: E402

IMAGES = Path(__file__).resolve().parent / "images"
IMAGES.mkdir(parents=True, exist_ok=True)
EQ_KW = dict(method="newton", tol=1e-10, max_iter=8000, min_concentration=1e-20)


def _ka(acid: str, base: str, pka: float) -> Reaction:
    return Reaction.from_string(f"{acid} > H+ & {base}", K=10 ** (-pka))


def selenium_environment(c_tot: float = 1.0) -> Enviroment:
    water = Compound("H2O", phase_point_list=[{"phase": "l", "temperature": 298}])
    kw = Reaction.from_string(
        f"{water.token} > H+ & OH-",
        concentrations={"H2O": XS(0.0)},
        K=1e-14,
    )
    return Enviroment(
        kw,
        _ka("HSeO4-", "SeO4-2", 1.92),
        _ka("H2SeO3", "HSeO3-", 2.62),
        _ka("HSeO3-", "SeO3-2", 7.19),
        _ka("H2Se", "HSe-", 3.89),
        HalfReaction.from_string(
            "HSeO4- & 3_H+ & 2_@e = H2SeO3 & H2O.l",
            E0=1.15,
            name="HSeO4-/H2SeO3",
        ),
        HalfReaction.from_string(
            "H2SeO3 & 4_H+ & 4_@e = Se.s & 3_H2O.l",
            E0=0.74,
            name="H2SeO3/Se",
        ),
        HalfReaction.from_string(
            "Se.s & 2_H+ & 2_@e = H2Se",
            E0=-0.11,
            name="Se/H2Se",
        ),
        concentrations={"HSeO4-": c_tot},
        buffer=["H+"],
    )


def write_pourbaix():
    diagram = Pourbaix(
        selenium_environment(),
        pH_min=0.0,
        pH_max=10.0,
        pH_steps=80,
        Eh_min=-1.0,
        Eh_max=1.4,
        Eh_steps=80,
        speciation_method="model",
    ).run()
    diagram.plot(
        plot_style="filled",
        save=str(IMAGES / "selenium_filled.png"),
        show=False,
    )
    diagram.plot(
        plot_style="labeled",
        save=str(IMAGES / "selenium_labeled.png"),
        show=False,
    )
    diagram.plot_boundaries(
        save=str(IMAGES / "selenium_boundaries.png"),
        show=False,
    )
    diagram.plot_predominance(
        save=str(IMAGES / "selenium_predominance.png"),
        show=False,
    )
    diagram.plot(
        plot_style="filled",
        show_frame_intersections=True,
        save=str(IMAGES / "selenium_frame.png"),
        show=False,
    )
    diagram.plot(
        plot_style="labeled",
        boundary_mode="all",
        save=str(IMAGES / "selenium_all_boundaries.png"),
        show=False,
    )
    plt.close("all")


def write_titration():
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
    curve = Titration(agf, titrant, volumes=volumes_ul * 1e-6).run(**EQ_KW)
    ksp = 4e-12
    qsp = curve.species("Ag+") ** 2 * curve.species("T-2")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(volumes_ul, qsp, color="#26547c", label=r"$Q_\mathrm{sp}=[\mathrm{Ag}^+]^2[\mathrm{T}^{2-}]$")
    ax.axhline(ksp, color="#ef476f", linestyle="--", label=r"$K_\mathrm{sp}(\mathrm{Ag}_2\mathrm{T})=4\times10^{-12}$")
    ax.set_xlabel(r"Titrant volume added ($\mu$L)")
    ax.set_ylabel(r"$Q_\mathrm{sp}$")
    ax.set_yscale("log")
    ax.legend()
    fig.tight_layout()
    fig.savefig(IMAGES / "ag2t_titration_qsp.png", dpi=150)
    plt.close("all")


def write_env16():
    env = build_env16()
    result = env.equilibrium(**EQ_KW, return_details=True)
    skip = {"H2O", "CaF2"}
    labels = []
    values = []
    for formula, value in result.concentrations_dict.items():
        if formula in skip:
            continue
        labels.append(formula)
        values.append(max(float(value), 1e-16))
    fig, ax = plt.subplots(figsize=(7, 4.8))
    y = np.arange(len(labels))
    ax.barh(y, values, color="#26547c")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xscale("log")
    ax.set_xlabel("Equilibrium concentration (mol/L)")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(IMAGES / "env16_speciation.png", dpi=150)
    plt.close("all")


def write_jungle():
    env = Enviroment(
        Reaction.from_string("R & T > 2_R", kf=0.01, kb=0.0, K=1e12),
        Reaction.from_string("R > D", kf=0.5, kb=0.0, K=1e12),
        Reaction.from_string("R & W > 2_W", kf=0.01, kb=0.0, K=1e12),
        Reaction.from_string("W > inert", kf=0.5, kb=0.0, K=1e12),
        concentrations={"R": 70.0, "T": 100.0, "W": 20.0, "D": 0.0, "inert": 0.0},
    )
    env.kinetics(
        time=40.0,
        accuracy=0.05,
        plot="save",
        directory=str(IMAGES / "jungle_model.png"),
        colors=["#2a9d8f", "#e9c46a", "#6d6875", "#e76f51", "#264653"],
    )


if __name__ == "__main__":
    write_pourbaix()
    write_titration()
    write_env16()
    write_jungle()
    print(f"Wrote figures to {IMAGES}")
