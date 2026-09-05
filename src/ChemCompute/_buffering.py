"""Solver-side environment buffering (constant species concentrations)."""


def clamp_buffered_concentrations(concentrations, env):
    """Reset buffered species to their fixed targets in-place."""
    targets = getattr(env, "_buffer_targets", None)
    if not targets:
        return concentrations
    labels = env.compound_labels
    for formula, target in targets.items():
        concentrations[labels.index(formula)] = target
    return concentrations


def merge_buffer_specs(*envs):
    """Union buffered species formulas from environments for combine."""
    merged = {}
    for env in envs:
        spec = getattr(env, "_buffer_spec", None) or {}
        for formula in spec:
            merged[formula] = None
    return merged
