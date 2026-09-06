"""Half-reactions, Nernst electrode potential, and redox equilibrium coupling."""

from __future__ import annotations

import copy
import math
import re
import warnings
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

R_GAS = 8.314462618  # J/(mol·K)
F_FARADAY = 96485.33212  # C/mol
LN10 = math.log(10.0)

_ELECTRON_ERROR = "Electrons belong in HalfReaction strings as @e, not in Reaction."


def _reject_electron_formula(formula: str) -> None:
    normalized = str(formula).strip().replace(" ", "")
    if normalized in {"e", "e-", "e+"}:
        raise ValueError(_ELECTRON_ERROR)


def _normalize_compound(entry):
    from ._formula import compound_from_species_token

    compound = entry["compound"]
    if isinstance(compound, str):
        _reject_electron_formula(compound)
        entry = dict(entry)
        entry["compound"] = compound_from_species_token(compound)
    elif hasattr(compound, "formula"):
        _reject_electron_formula(compound.formula)
    stoich = float(entry.get("stoichiometric_coefficient", 1))
    entry.setdefault("stoichiometric_coefficient", stoich)
    entry.setdefault("rate_dependency", stoich)
    return entry


def _stoichiometry_key(species_list: list[dict]) -> tuple:
    items = []
    for entry in species_list:
        compound = entry["compound"]
        formula = compound.formula if hasattr(compound, "formula") else str(compound)
        coeff = float(entry.get("stoichiometric_coefficient", 1))
        items.append((formula, coeff))
    return tuple(sorted(items))


def compute_pH(env, concentrations: Optional[np.ndarray] = None) -> float:
    """Return pH from env H+ concentration (activity-corrected when model active)."""
    conc = np.asarray(
        concentrations if concentrations is not None else env.concentrations,
        dtype=float,
    )
    h_index = None
    for j, compound in enumerate(env.compounds):
        if compound.formula == "H+":
            h_index = j
            break
    if h_index is None:
        return 7.0
    h_conc = conc[h_index]
    if getattr(env, "activity_model", None) is not None:
        gammas = env.activity_model.gamma_array(env, np.maximum(conc, 0.0), env.T)
        h_conc = gammas[h_index] * h_conc
    return -math.log10(max(h_conc, 1e-300))


def redox_lnK(hr, Eh: float, pH: float, T: float) -> float:
    """Natural log of equilibrium K for a half-reaction at electrode potential Eh."""
    return math.log(max(hr.K_at(Eh, pH, T=T), 1e-300))


def d_redox_lnK_dEh(hr, T: float) -> float:
    """Derivative d(ln K)/dEh for Nernst K_at."""
    n = hr.n_electrons
    return n * F_FARADAY / (R_GAS * T)


def apply_electrode_potential(env) -> None:
    """Mode A: overwrite K on linked half-reactions when env.electrode_Eh is set."""
    eh = getattr(env, "electrode_Eh", None)
    half_reactions = getattr(env, "half_reactions", None) or []
    if eh is None or not half_reactions:
        return
    pH = compute_pH(env)
    for hr in half_reactions:
        idx = hr._reaction_index
        if idx is None or idx >= len(env.reactions):
            continue
        rxn = env.reactions[idx]
        if not hasattr(rxn, "_K_ref"):
            rxn._K_ref = rxn.K
        rxn.K = hr.K_at(eh, pH, T=env.T)


def _tokenize_half_side(side: str) -> list[str]:
    side = side.replace(" ", "")
    side = re.sub(r"\+(?=@e)", "", side)
    electron = r"(?:\d+(?:\.\d+)?_)?@e"
    live = r"(?:\d+(?:\.\d+)?_)?@c\d+(?:_-?\d+(?:\.\d+)?)?"
    ion_suffix = r"(?:[+-]\d*|\+\d*|-\d*)?"
    term = (
        rf"\d+(?:\.\d+)?_[A-Za-z0-9.\-()\[\]]+{ion_suffix}(?:\.(?:s|l|g|aq))?"
        rf"|[A-Za-z0-9.\-()\[\]]+{ion_suffix}(?:\.(?:s|l|g|aq))?"
    )
    pattern = rf"(?:{electron}|{live}|{term})"
    tokens = re.findall(pattern, side)
    if not tokens:
        raise ValueError(f"Could not parse half-reaction side: {side!r}")
    return tokens


def _split_side_sections(side: str) -> list[str]:
    side = side.replace(" ", "")
    chunks = [chunk for chunk in side.split("&") if chunk]
    tokens: list[str] = []
    for chunk in chunks:
        tokens.extend(_tokenize_half_side(chunk))
    return tokens


def _parse_side_tokens(side: str) -> tuple[list[dict], float]:
    """Parse one side of a half-reaction string; return species list and electron count."""
    sections = _split_side_sections(side)
    species: list[dict] = []
    n_electrons = 0.0
    for section in sections:
        if section == "@e":
            n_electrons += 1.0
            continue
        electron_match = re.match(r"^(\d+(?:\.\d+)?)_@e$", section)
        if electron_match:
            n_electrons += float(electron_match.group(1))
            continue
        if section in {"e", "e-", "e+"}:
            raise ValueError(
                "Use @e as the electron token in half-reaction strings, not bare e or e-."
            )
        if "@" in section:
            from ._interpolation import is_live_compound_section

            if not is_live_compound_section(section):
                raise ValueError(f"Invalid electron token {section!r}; use @e or n_@e.")

        species.append(_parse_section(section))
    return species, n_electrons


def _parse_section(section: str) -> dict:
    _species = r"(?:@c\d+|[A-Za-z0-9+.\-()\[\]]+)"
    acceptable = re.compile(
        rf"^(\d+(?:\.\d+)?_{_species}_-?\d+(?:\.\d+)?|"
        rf"\d+(?:\.\d+)?_{_species}|"
        rf"{_species}_-?\d+(?:\.\d+)?|"
        rf"{_species})(\.s|\.g|\.l|\.aq)?$"
    )
    if not acceptable.match(section):
        raise ValueError(f"Invalid half-reaction term: {section!r}")

    from ._formula import split_phase_suffix

    core, phase = split_phase_suffix(section)
    phase_suffix = f".{phase}" if phase else ""

    parts = core.split("_")
    if len(parts) == 3:
        info = {
            "stoichiometric_coefficient": float(parts[0]),
            "compound": parts[1] + phase_suffix,
            "rate_dependency": float(parts[2]),
        }
    elif len(parts) == 2:
        if re.match(r"^\d+(?:\.\d+)?$", parts[0]):
            stoich = float(parts[0])
            info = {
                "stoichiometric_coefficient": stoich,
                "compound": parts[1] + phase_suffix,
                "rate_dependency": stoich,
            }
        else:
            info = {
                "stoichiometric_coefficient": 1,
                "compound": parts[0] + phase_suffix,
                "rate_dependency": float(parts[1]),
            }
    else:
        info = {
            "stoichiometric_coefficient": 1,
            "compound": parts[0] + phase_suffix,
            "rate_dependency": 1,
        }
    return info


def _assign_compounds(species_list: list[dict], T: float) -> None:
    from ._general import Compound
    from ._interpolation import compound_from_parsed_name

    for entry in species_list:
        name = entry["compound"]
        if isinstance(name, Compound):
            continue
        entry["compound"] = compound_from_parsed_name(name, T=T)


@dataclass
class BoundaryLine:
    """Analytical Pourbaix boundary line for one half-reaction."""

    half_reaction: "HalfReaction"
    pH_values: np.ndarray
    Eh_values: np.ndarray
    label: str = ""


@dataclass
class HalfReaction:
    """
    Electrochemical half-reaction with Nernst electrode potential coupling.

    Species lists mirror :class:`Reaction` reactants/products. Use ``@e`` in string
    parsers as the electron token (never ``e`` or ``e-``).
    """

    oxidized: list[dict]
    reduced: list[dict]
    oxidized_concentration: list[float]
    reduced_concentration: list[float]
    E0_SHE: float
    n_electrons: float = 1.0
    T: float = 298.0
    name: str = ""
    compounds: list[dict] = field(default_factory=list)
    _reaction_index: Optional[int] = None

    def __post_init__(self):
        self.oxidized = [_normalize_compound(dict(e)) for e in self.oxidized]
        self.reduced = [_normalize_compound(dict(e)) for e in self.reduced]
        self._build_compounds_list()

    def _build_compounds_list(self) -> None:
        self.compounds = []
        for entry, conc in zip(self.oxidized, self.oxidized_concentration):
            item = dict(entry)
            item["concentration"] = conc
            item["type"] = "oxidized"
            self.compounds.append(item)
        for entry, conc in zip(self.reduced, self.reduced_concentration):
            item = dict(entry)
            item["concentration"] = conc
            item["type"] = "reduced"
            self.compounds.append(item)

    @classmethod
    def from_string(
        cls,
        reaction_str: str,
        *,
        concentrations: Optional[list[float]] = None,
        E0: float = 0.0,
        T: float = 298,
        name: str = "",
    ) -> "HalfReaction":
        if "=" not in reaction_str:
            raise ValueError("Half-reaction strings must use '=' as the separator.")
        left, right = reaction_str.replace(" ", "").split("=", 1)
        oxidized, n_left = _parse_side_tokens(left)
        reduced, n_right = _parse_side_tokens(right)
        n_electrons = n_left + n_right
        if n_electrons <= 0:
            raise ValueError("Half-reaction must include @e on the oxidized (left) side.")
        if n_right > 0:
            raise ValueError("@e must appear only on the oxidized (left) side of '='.")

        _assign_compounds(oxidized, T)
        _assign_compounds(reduced, T)

        n_species = len(oxidized) + len(reduced)
        if concentrations is None:
            concentrations = [0.0] * n_species
        if len(concentrations) != n_species:
            raise ValueError(
                f"Expected {n_species} concentrations (excluding @e), got {len(concentrations)}."
            )
        ox_conc = concentrations[: len(oxidized)]
        red_conc = concentrations[len(oxidized) :]
        return cls(
            oxidized=oxidized,
            reduced=reduced,
            oxidized_concentration=list(ox_conc),
            reduced_concentration=list(red_conc),
            E0_SHE=float(E0),
            n_electrons=float(n_electrons),
            T=float(T),
            name=name or reaction_str,
        )

    def net_h_plus_stoichiometry(self) -> float:
        """Signed m in Ox + m H+ ⇌ Red (positive when H+ is consumed)."""
        m = 0.0
        for entry in self.oxidized:
            if entry["compound"].formula == "H+":
                m += float(entry["stoichiometric_coefficient"])
        for entry in self.reduced:
            if entry["compound"].formula == "H+":
                m -= float(entry["stoichiometric_coefficient"])
        return m

    def reduction_stoichiometry(self) -> tuple[list[dict], list[dict]]:
        """Mass-action reactants/products without electrons."""
        reactants = copy.deepcopy(self.oxidized)
        products = copy.deepcopy(self.reduced)
        return reactants, products

    def E_at_pH(self, pH: float, T: Optional[float] = None) -> float:
        """Standard Nernst line at unit activity vs pH."""
        T = self.T if T is None else T
        n = self.n_electrons
        if n <= 0:
            raise ValueError("n_electrons must be positive.")
        m = self.net_h_plus_stoichiometry()
        return self.E0_SHE - (R_GAS * T / (n * F_FARADAY)) * LN10 * m * pH

    def K_at(self, Eh: float, pH: float, T: Optional[float] = None) -> float:
        """Equilibrium constant at electrode potential Eh (V vs SHE)."""
        T = self.T if T is None else T
        n = self.n_electrons
        e_at = self.E_at_pH(pH, T=T)
        return math.exp(n * F_FARADAY * (Eh - e_at) / (R_GAS * T))

    def _concentration_vector(self, env=None, concentrations: Optional[np.ndarray] = None) -> np.ndarray:
        if env is not None and self._reaction_index is not None:
            return np.asarray(env.concentrations, dtype=float)
        values = self.oxidized_concentration + self.reduced_concentration
        return np.asarray(values, dtype=float)

    def _build_local_env(self, concentrations: Optional[np.ndarray] = None):
        """Standalone Reaction/env for Q when not attached."""
        from ._general import Enviroment, Reaction

        reactants, products = self.reduction_stoichiometry()
        rxn = Reaction(
            reactants,
            products,
            list(self.oxidized_concentration),
            list(self.reduced_concentration),
            K=1.0,
            T=self.T,
        )
        return Enviroment(rxn, T=self.T)

    def reaction_quotient(
        self,
        env=None,
        concentrations: Optional[np.ndarray] = None,
    ) -> float:
        """Mass-action quotient Q using env activities or stored concentrations."""
        from ._equilibrium import _build_context, _compute_lnQ, _safe_concentrations

        if env is None:
            env = self._build_local_env(concentrations)
        elif concentrations is not None:
            conc = np.asarray(concentrations, dtype=float)
            for j, val in enumerate(conc):
                if j < len(env.compounds_concentration):
                    env.compounds_concentration[j]["concentration"] = float(val)

        ctx = _build_context(env, min_concentration=1e-12)
        idx = self._reaction_index if self._reaction_index is not None else 0
        c = np.asarray(env.concentrations, dtype=float)
        c_safe = _safe_concentrations(c, 1e-12)
        lnQ_all = _compute_lnQ(ctx, c_safe, c)
        return float(math.exp(lnQ_all[idx]))

    def E_at(
        self,
        env=None,
        concentrations: Optional[np.ndarray] = None,
        T: Optional[float] = None,
    ) -> float:
        """Nernst electrode potential (V vs SHE) from current Q."""
        T = self.T if T is None else T
        n = self.n_electrons
        q = self.reaction_quotient(env=env, concentrations=concentrations)
        q = max(q, 1e-300)
        return self.E0_SHE - (R_GAS * T / (n * F_FARADAY)) * math.log(q)

    def boundary_line(self, pH_values: Sequence[float]) -> BoundaryLine:
        pH_arr = np.asarray(pH_values, dtype=float)
        eh = np.array([self.E_at_pH(float(p)) for p in pH_arr])
        return BoundaryLine(
            half_reaction=self,
            pH_values=pH_arr,
            Eh_values=eh,
            label=self.name or "half-reaction",
        )

    def _linked_reaction(self):
        from ._general import Reaction

        reactants, products = self.reduction_stoichiometry()
        return Reaction(
            reactants,
            products,
            list(self.oxidized_concentration),
            list(self.reduced_concentration),
            K=1.0,
            T=self.T,
        )

    def attach(self, env) -> int:
        """Link to an existing Reaction with matching stoichiometry; return its index."""
        target_r, target_p = self.reduction_stoichiometry()
        key_r = _stoichiometry_key(target_r)
        key_p = _stoichiometry_key(target_p)
        for i, rxn in enumerate(env.reactions):
            if _stoichiometry_key(rxn.reactants) == key_r and _stoichiometry_key(rxn.products) == key_p:
                self._reaction_index = i
                return i
        raise ValueError("No matching Reaction stoichiometry found in environment.")

    def attach_or_create(self, env) -> int:
        """Ensure a linked Reaction exists in env; return its index."""
        try:
            return self.attach(env)
        except ValueError:
            pass
        rxn = self._linked_reaction()
        rxn._adjust_thermodynamics = env.adjust_thermodynamics
        rxn.T = env.T
        env.reactions.append(rxn)
        self._reaction_index = len(env.reactions) - 1
        return self._reaction_index

    def stoichiometry_signature(self) -> tuple:
        reactants, products = self.reduction_stoichiometry()
        return (_stoichiometry_key(reactants), _stoichiometry_key(products))


def register_half_reactions(env, half_reactions: Sequence[HalfReaction]) -> None:
    """Register half-reactions on env with duplicate validation."""
    seen: dict[tuple, HalfReaction] = {}
    for existing in getattr(env, "half_reactions", None) or []:
        seen[existing.stoichiometry_signature()] = existing
    registered: list[HalfReaction] = list(getattr(env, "half_reactions", None) or [])
    for hr in half_reactions:
        sig = hr.stoichiometry_signature()
        if sig in seen and seen[sig].E0_SHE != hr.E0_SHE:
            raise ValueError(
                "Duplicate half-reaction stoichiometry with different E0 values."
            )
        seen[sig] = hr
        if hr not in registered:
            registered.append(hr)
    env.half_reactions = registered


def redox_reaction_indices(env) -> list[int]:
    return [
        hr._reaction_index
        for hr in (getattr(env, "half_reactions", None) or [])
        if hr._reaction_index is not None
    ]


def coupled_eh_mode(env) -> bool:
    hrs = getattr(env, "half_reactions", None) or []
    return len(hrs) >= 2 and getattr(env, "electrode_Eh", None) is None


def update_redox_lnK(ctx, Eh: float) -> None:
    """Refresh lnK on redox rows of an equilibrium context."""
    env = ctx.env
    pH = compute_pH(env, ctx.c0 + ctx.S @ getattr(ctx, "_last_x", np.zeros(ctx.R)))
    for hr in getattr(env, "half_reactions", None) or []:
        idx = hr._reaction_index
        if idx is None:
            continue
        ctx.lnK[idx] = redox_lnK(hr, Eh, pH, env.T)


def initial_eh_guess(env) -> float:
    """Average E_at over registered half-reactions at current concentrations."""
    hrs = getattr(env, "half_reactions", None) or []
    if not hrs:
        return 0.0
    values = [hr.E_at(env) for hr in hrs]
    return float(np.mean(values))
