"""Epidemic population structure simulation package."""

from .contact_matrix import build_contact_matrix, preset_general_cluster
from .interventions import (
    ContactMatrixIntervention,
    InterventionPlan,
    compose_intervention_plan,
    intervene_bridge_links,
    intervene_within_groups,
)
from .model import (
    EpidemicModel,
    RK4Stepper,
    StatePoint,
    StateStepper,
    StructuredEpidemicModel,
    StructuredStatePoint,
)
from .parameters import DiseaseParameters, PopulationParameters, SimulationParameters
from .scenarios import Scenario, get_scenario, list_scenarios
from .notebook_helpers import (
    draw_group_bridge_map,
    effective_matrix,
    get_contact_matrix,
    get_intervention_plan,
    intervention_targets,
    plot_matrix_transformation,
    plot_population_infected_comparison,
    run_scenario,
    summarize_two_population_runs,
)

__all__ = [
    "build_contact_matrix",
    "ContactMatrixIntervention",
    "compose_intervention_plan",
    "DiseaseParameters",
    "draw_group_bridge_map",
    "effective_matrix",
    "EpidemicModel",
    "get_contact_matrix",
    "get_intervention_plan",
    "get_scenario",
    "intervention_targets",
    "InterventionPlan",
    "intervene_bridge_links",
    "intervene_within_groups",
    "list_scenarios",
    "plot_matrix_transformation",
    "plot_population_infected_comparison",
    "PopulationParameters",
    "preset_general_cluster",
    "RK4Stepper",
    "run_scenario",
    "Scenario",
    "SimulationParameters",
    "StatePoint",
    "StateStepper",
    "StructuredEpidemicModel",
    "StructuredStatePoint",
    "summarize_two_population_runs",
]
__version__ = "0.1.0"
