"""Reusable scenario definitions for epidemiological simulations.

Each scenario encapsulates a complete setup: populations, contact matrix,
disease parameters, simulation configuration, and optional interventions.

Scenarios are used by notebooks to avoid repeated parameter construction
and to maintain canonical baseline + intervention variants for comparison.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd

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
    """A complete simulation scenario: populations, contact matrix, and interventions.

    Common usage patterns:

    - Build variants from an existing scenario with ``Scenario.builder(base=...)``.
    - Read one transition using tuple indexing, e.g. ``scenario["general", "cluster"]``.
    - Manually tweak one transition using assignment,
        e.g. ``scenario["general", "cluster"] = 0.0001``.
    """

    name: str
    description: str
    populations: list[PopulationParameters]
    contact_matrix: list[list[float]]
    population_names: list[str]
    simulation: SimulationParameters
    intervention_plan: Optional[InterventionPlan] = None

    def transition(self, recipient_group: str | int, source_group: str | int) -> float:
        """Return one contact-matrix entry by name or index."""

        recipient_index = _resolve_group_ref(self.population_names, recipient_group)
        source_index = _resolve_group_ref(self.population_names, source_group)
        return float(self.contact_matrix[recipient_index][source_index])

    def set_transition(
        self,
        recipient_group: str | int,
        source_group: str | int,
        value: float,
    ) -> None:
        """Set one contact-matrix entry by name or index."""

        recipient_index = _resolve_group_ref(self.population_names, recipient_group)
        source_index = _resolve_group_ref(self.population_names, source_group)
        self.contact_matrix[recipient_index][source_index] = float(value)

    def __getitem__(self, key: tuple[str | int, str | int]) -> float:
        """Convenience lookup: scenario[recipient, source]."""

        if not isinstance(key, tuple) or len(key) != 2:
            raise TypeError("Scenario transition lookup expects scenario[recipient, source]")
        recipient_group, source_group = key
        return self.transition(recipient_group, source_group)

    def __setitem__(self, key: tuple[str | int, str | int], value: float) -> None:
        """Convenience assignment: scenario[recipient, source] = value."""

        if not isinstance(key, tuple) or len(key) != 2:
            raise TypeError("Scenario transition assignment expects scenario[recipient, source] = value")
        recipient_group, source_group = key
        self.set_transition(recipient_group, source_group, value)

    @classmethod
    def builder(
        cls,
        base: Scenario | str | None = None,
        *,
        name: str | None = None,
        description: str | None = None,
        simulation: SimulationParameters | None = None,
        intervention_plan: InterventionPlan | None = None,
    ) -> ScenarioBuilder:
        """Create a builder seeded from a base scenario.

        If ``base`` is omitted, the baseline scenario is used as the default template.
        """

        return ScenarioBuilder(
            base=base,
            name=name,
            description=description,
            simulation=simulation,
            intervention_plan=intervention_plan,
        )
    
    def contact_matrix_df(self) -> pd.DataFrame:
        """Return the contact matrix as a pandas DataFrame with population names as index/columns."""
        return pd.DataFrame(self.contact_matrix, index=self._population_names, columns=self._population_names)


class ScenarioBuilder:
    """Builder for assembling a scenario incrementally.

    The builder is centered on a Scenario so callers can start from a base
    template, override selected defaults, and define populations and contact
    matrix entries one at a time.
    """

    def __init__(
        self,
        base: Scenario | str | None = None,
        *,
        name: str | None = None,
        description: str | None = None,
        simulation: SimulationParameters | None = None,
        intervention_plan: InterventionPlan | None = None,
    ) -> None:
        base_scenario = _resolve_base_scenario(base)

        self._name = name or base_scenario.name
        self._description = description or base_scenario.description
        self._simulation = simulation or base_scenario.simulation
        self._intervention_plan = (
            base_scenario.intervention_plan
            if intervention_plan is None
            else intervention_plan
        )
        self._populations = list(base_scenario.populations)
        self._population_names = list(base_scenario.population_names)
        self._population_index = {
            population_name: index
            for index, population_name in enumerate(self._population_names)
        }
        self._contact_matrix = [row[:] for row in base_scenario.contact_matrix]
        self._default_disease = _infer_default_disease(base_scenario)

    def name(self, name: str) -> ScenarioBuilder:
        """Override the scenario name."""

        self._name = name
        return self

    def description(self, description: str) -> ScenarioBuilder:
        """Override the scenario description."""

        self._description = description
        return self

    def disease(self, disease: DiseaseParameters) -> ScenarioBuilder:
        """Set the default disease parameters used by later population entries."""

        self._default_disease = disease
        return self

    def simulation(self, simulation: SimulationParameters) -> ScenarioBuilder:
        """Override the simulation controls."""

        self._simulation = simulation
        return self

    def intervention_plan(self, intervention_plan: InterventionPlan | None) -> ScenarioBuilder:
        """Attach or clear the scenario intervention calendar."""

        self._intervention_plan = intervention_plan
        return self

    def population(
        self,
        population: PopulationParameters | str,
        *,
        size: int | None = None,
        beta: float | None = None,
        initial_susceptible: int | None = None,
        initial_exposed: int | None = None,
        initial_infected: int | None = None,
        initial_recovered: int | None = None,
        disease: DiseaseParameters | None = None,
    ) -> ScenarioBuilder:
        """Register or replace a population definition."""

        if isinstance(population, PopulationParameters):
            population_parameters = population
        else:
            existing = self._population_by_name(population)
            if existing is None and size is None:
                raise ValueError("size must be provided when registering a new population")

            population_size = size if size is not None else existing.size
            population_beta = (
                existing.beta if existing is not None and beta is None else (0.0 if beta is None else beta)
            )
            population_susceptible = (
                existing.initial_susceptible
                if existing is not None and initial_susceptible is None
                else initial_susceptible
            )
            population_exposed = (
                existing.initial_exposed
                if existing is not None and initial_exposed is None
                else (0 if initial_exposed is None else initial_exposed)
            )
            population_infected = (
                existing.initial_infected
                if existing is not None and initial_infected is None
                else (0 if initial_infected is None else initial_infected)
            )
            population_recovered = (
                existing.initial_recovered
                if existing is not None and initial_recovered is None
                else (0 if initial_recovered is None else initial_recovered)
            )
            population_disease = (
                disease
                or (existing.disease if existing is not None else None)
                or self._default_disease
                or DiseaseParameters()
            )

            population_parameters = PopulationParameters(
                name=population,
                size=population_size,
                beta=population_beta,
                initial_susceptible=population_susceptible,
                initial_exposed=population_exposed,
                initial_infected=population_infected,
                initial_recovered=population_recovered,
                disease=population_disease,
            )

        self._register_population(population_parameters)
        return self

    def transition(
        self,
        recipient_group: str | int,
        source_group: str | int,
        value: float,
    ) -> ScenarioBuilder:
        """Set a single contact-matrix entry.

        Group identifiers may be names or zero-based indices, but the group must
        already be registered in the scenario.
        """

        recipient_index = self._resolve_group(recipient_group)
        source_index = self._resolve_group(source_group)
        self._contact_matrix[recipient_index][source_index] = float(value)
        return self

    def within(self, group: str | int, value: float) -> ScenarioBuilder:
        """Set a within-group transition on the matrix diagonal."""

        return self.transition(group, group, value)

    def bridge(self, recipient_group: str | int, source_group: str | int, value: float) -> ScenarioBuilder:
        """Set an off-diagonal bridge transition."""

        if recipient_group == source_group:
            raise ValueError("bridge transitions must be off-diagonal")
        return self.transition(recipient_group, source_group, value)

    def build(self) -> Scenario:
        """Finalize the builder as an immutable Scenario."""

        return Scenario(
            name=self._name,
            description=self._description,
            populations=list(self._populations),
            contact_matrix=[row[:] for row in self._contact_matrix],
            population_names=list(self._population_names),
            simulation=self._simulation,
            intervention_plan=self._intervention_plan,
        )

    def _register_population(self, population: PopulationParameters) -> None:
        name = population.name
        if name in self._population_index:
            index = self._population_index[name]
            self._populations[index] = population
            return

        self._population_index[name] = len(self._populations)
        self._population_names.append(name)
        self._populations.append(population)

        for row in self._contact_matrix:
            row.append(0.0)
        self._contact_matrix.append([0.0 for _ in range(len(self._population_names))])

    def _population_by_name(self, name: str) -> PopulationParameters | None:
        index = self._population_index.get(name)
        if index is None:
            return None
        return self._populations[index]

    def _resolve_group(self, group: str | int) -> int:
        return _resolve_group_ref(self._population_names, group)


def _resolve_base_scenario(base: Scenario | str | None) -> Scenario:
    if base is None:
        return get_scenario("baseline")
    if isinstance(base, str):
        return get_scenario(base)
    return base


def _infer_default_disease(base: Scenario) -> DiseaseParameters:
    if not base.populations:
        return DiseaseParameters()
    return base.populations[0].disease


def _resolve_group_ref(population_names: list[str], group: str | int) -> int:
    if isinstance(group, int):
        if group < 0 or group >= len(population_names):
            raise ValueError(f"Unknown population index: {group}")
        return group

    if group not in population_names:
        valid = ", ".join(population_names)
        raise ValueError(f"Unknown population '{group}'. Registered populations: {valid}")
    return population_names.index(group)


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
