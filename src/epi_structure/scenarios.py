"""Reusable scenario definitions for epidemiological simulations.

Each scenario encapsulates a complete setup: populations, contact matrix,
disease parameters, simulation configuration, and optional interventions.

Scenarios are used by notebooks to avoid repeated parameter construction
and to maintain canonical baseline + intervention variants for comparison.
"""

from dataclasses import dataclass
from typing import Optional

from .contact_matrix import preset_general_cluster
from .interventions import (
    InterventionPlan,
    compose_intervention_plan,
    intervene_bridge_links,
    intervene_within_groups,
)
from .parameters import (
    DiseaseParameters,
    PopulationParameters,
    SimulationParameters,
)


@dataclass
class Scenario:
    """A complete simulation scenario: populations, contact matrix, and interventions."""

    name: str
    description: str
    populations: list[PopulationParameters]
    contact_matrix: list[list[float]]
    population_names: list[str]
    simulation: SimulationParameters
    intervention_plan: Optional[InterventionPlan] = None


def _baseline_scenario() -> Scenario:
    """Baseline (no intervention) two-population SEIRS scenario.

    Setup:
    - General population: 10,000 individuals, 1 initial infected
    - Cluster population: 500 individuals, 4 initial infected
    - Disease: infectious period 6 days, latent period 2.5 days
    - Duration: 60 days, output every 1 day
    - Contact matrix: within-group and bridging transmission
    """
    shared_disease = DiseaseParameters(
        infectious_period=6.0, latent_period=2.5
    )
    populations = [
        PopulationParameters(
            name="general",
            size=10_000,
            beta=0.0,
            initial_infected=1,
            disease=shared_disease,
        ),
        PopulationParameters(
            name="cluster",
            size=500,
            beta=0.0,
            initial_infected=4,
            disease=shared_disease,
        ),
    ]
    names, contact_matrix = preset_general_cluster()
    simulation = SimulationParameters(
        time_step=0.1, duration=60.0, output_stride=10
    )
    return Scenario(
        name="baseline",
        description="No intervention (natural progression)",
        populations=populations,
        contact_matrix=contact_matrix,
        population_names=names,
        simulation=simulation,
        intervention_plan=None,
    )


def _cluster_ventilation_scenario() -> Scenario:
    """Cluster ventilation intervention starting at t=20.

    Policy: Improve ventilation/masking in cluster settings → scale cluster
    within-group transmission by 40% (scale=0.4) beginning at t=20.
    """
    base = _baseline_scenario()
    intervention_plan = compose_intervention_plan(
        intervene_within_groups(["cluster"], 0.4, start_time=20.0)
    )
    return Scenario(
        name="cluster_ventilation",
        description="Ventilation/masking reduces cluster within-group transmission to 40% from t=20",
        populations=base.populations,
        contact_matrix=base.contact_matrix,
        population_names=base.population_names,
        simulation=base.simulation,
        intervention_plan=intervention_plan,
    )


def _bridge_tracing_scenario() -> Scenario:
    """Contact tracing on bridge links starting at t=20.

    Policy: Contact tracing + isolation of bridge cases → scale bidirectional
    bridge transmission between general and cluster by 20% (scale=0.2) from t=20.
    """
    base = _baseline_scenario()
    intervention_plan = compose_intervention_plan(
        intervene_bridge_links(
            [("general", "cluster")],
            0.2,
            start_time=20.0,
            symmetric=True,
        )
    )
    return Scenario(
        name="bridge_tracing",
        description="Contact tracing reduces bridge transmission to 20% from t=20",
        populations=base.populations,
        contact_matrix=base.contact_matrix,
        population_names=base.population_names,
        simulation=base.simulation,
        intervention_plan=intervention_plan,
    )


def _targeted_combo_scenario() -> Scenario:
    """Combined intervention: cluster ventilation + bridge tracing from t=20.

    Policy: Apply both cluster within-group reduction (40%) and bridge
    tracing (20%) simultaneously, starting at t=20.
    """
    base = _baseline_scenario()
    intervention_plan = compose_intervention_plan(
        intervene_within_groups(["cluster"], 0.4, start_time=20.0),
        intervene_bridge_links(
            [("general", "cluster")],
            0.2,
            start_time=20.0,
            symmetric=True,
        ),
    )
    return Scenario(
        name="targeted_combo",
        description="Ventilation (40%) + bridge tracing (20%) combined from t=20",
        populations=base.populations,
        contact_matrix=base.contact_matrix,
        population_names=base.population_names,
        simulation=base.simulation,
        intervention_plan=intervention_plan,
    )


# Registry of canonical scenarios
_SCENARIO_REGISTRY = {
    "baseline": _baseline_scenario,
    "cluster_ventilation": _cluster_ventilation_scenario,
    "bridge_tracing": _bridge_tracing_scenario,
    "targeted_combo": _targeted_combo_scenario,
}


def get_scenario(name: str) -> Scenario:
    """Retrieve a scenario by name.

    Args:
        name: Scenario identifier (e.g., "baseline", "cluster_ventilation")

    Returns:
        Scenario object with all parameters configured

    Raises:
        ValueError: If name is not in the scenario registry
    """
    if name not in _SCENARIO_REGISTRY:
        valid = ", ".join(sorted(_SCENARIO_REGISTRY.keys()))
        raise ValueError(
            f"Unknown scenario '{name}'. Valid options: {valid}"
        )
    return _SCENARIO_REGISTRY[name]()


def list_scenarios() -> list[str]:
    """Return list of available scenario names."""
    return sorted(_SCENARIO_REGISTRY.keys())
