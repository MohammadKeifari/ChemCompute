"""Volume-aware environment mixing with optional coefficients."""

from __future__ import annotations

import copy as copy_module


def _resolve_compound_key(key, T: float):
    from ._general import Compound

    if isinstance(key, Compound):
        return key
    return Compound(key, phase_point_list=[{"phase": "aq", "temperature": T}])


def _validate_coefficient(coefficient: float) -> float:
    coefficient = float(coefficient)
    if coefficient <= 0:
        raise ValueError("Mixing coefficient must be positive.")
    return coefficient


class ScaledEnviroment:
    """Wrapper for coefficient-weighted environment mixing (e.g. 0.5 * envA)."""

    def __init__(self, coefficient: float, env):
        from ._general import Enviroment

        self.coefficient = _validate_coefficient(coefficient)
        if not isinstance(env, Enviroment):
            raise TypeError("ScaledEnviroment requires an Enviroment instance.")
        self.env = env

    def __add__(self, other):
        from ._general import Enviroment

        if isinstance(other, ScaledEnviroment):
            return Enviroment.combine(self, other)
        if isinstance(other, Enviroment):
            return Enviroment.combine(self, (1.0, other))
        return NotImplemented

    def __radd__(self, other):
        if other == 0:
            return self
        return NotImplemented


def _normalize_combine_terms(*terms):
    from ._general import Enviroment

    normalized = []
    for term in terms:
        if isinstance(term, ScaledEnviroment):
            normalized.append((term.coefficient, term.env))
        elif isinstance(term, Enviroment):
            normalized.append((1.0, term))
        elif isinstance(term, tuple) and len(term) == 2:
            coeff, env = term
            if not isinstance(env, Enviroment):
                raise TypeError("Combine tuple must be (coefficient, Enviroment).")
            normalized.append((_validate_coefficient(coeff), env))
        else:
            raise TypeError(
                "Each combine term must be Enviroment, ScaledEnviroment, or (coeff, Enviroment)."
            )
    if not normalized:
        raise ValueError("At least one environment term is required for combine.")
    return normalized


def combine_environments(*terms):
    """Merge environments with coefficient-weighted volume mixing."""
    from ._general import Enviroment, Reaction

    normalized = _normalize_combine_terms(*terms)

    reference = normalized[0][1]
    T = reference.T
    adjust_thermodynamics = reference.adjust_thermodynamics
    activity_model = reference.activity_model

    for _, env in normalized[1:]:
        if env.T != T:
            raise ValueError("All environments must have the same temperature to combine.")
        if env.adjust_thermodynamics != adjust_thermodynamics:
            raise ValueError("All environments must have the same adjust_thermodynamics setting.")
        if env.activity_model is not None and activity_model is not None:
            if getattr(env.activity_model, "model", None) != getattr(activity_model, "model", None):
                raise ValueError("Conflicting activity models in combine.")
        if activity_model is None and env.activity_model is not None:
            activity_model = env.activity_model

    total_effective_volume = sum(coeff * env.volume for coeff, env in normalized)

    compound_objects = {}
    mole_totals = {}
    excess_flags = {}
    merged_reactions = []
    charge_map = {}

    for coeff, env in normalized:
        merged_reactions.extend(copy_module.deepcopy(env.reactions))
        charge_map.update(env.charge_map)
        effective_volume = coeff * env.volume
        env_excess = getattr(env, "excess_dict", {})
        for formula, concentration in env.concentrations_dict.items():
            incoming = env.compounds[env.compound_labels.index(formula)]
            if formula not in compound_objects:
                compound_objects[formula] = incoming
            elif getattr(incoming, "spectrum", None) is not None:
                compound_objects[formula] = incoming
            mole_totals[formula] = mole_totals.get(formula, 0.0) + concentration * effective_volume
            excess_flags[formula] = excess_flags.get(formula, False) or env_excess.get(formula, False)

    combined = Enviroment.__new__(Enviroment)
    combined.reactions = merged_reactions
    combined.adjust_thermodynamics = adjust_thermodynamics
    combined.charge_map = charge_map
    combined._activity_model = activity_model
    combined._T = T
    combined.volume = total_effective_volume
    combined.compounds = []
    combined.compounds_concentration = []
    combined._last_equilibrium_result = None

    for formula in sorted(compound_objects.keys()):
        combined.compounds.append(compound_objects[formula])
        combined.compounds_concentration.append(
            {
                "compound": compound_objects[formula],
                "concentration": mole_totals[formula] / total_effective_volume,
                "excess": excess_flags.get(formula, False),
            }
        )

    from ._buffering import merge_buffer_specs

    combined._buffer_spec = merge_buffer_specs(*(env for _, env in normalized))
    combined._resolve_buffer_targets()

    for reaction in combined.reactions:
        reaction._adjust_thermodynamics = adjust_thermodynamics
        reaction.T = T

    return combined


def add_compounds_to_environment(
    env,
    concentrations: dict,
    *,
    volume: float = 1.0,
    coefficient: float = 1.0,
):
    """Pour a concentration slug into an environment with optional volume coefficient."""
    from ._general import Enviroment

    if volume <= 0:
        raise ValueError("Slug volume must be positive.")
    slug = Enviroment.from_compounds(
        concentrations,
        T=env.T,
        volume=volume,
        adjust_thermodynamics=env.adjust_thermodynamics,
        activity_model=env.activity_model,
    )
    return combine_environments((coefficient, slug), (1.0, env))
