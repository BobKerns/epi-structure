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

__all__ = [
    "build_contact_matrix",
    "ContactMatrixIntervention",
    "compose_intervention_plan",
    "DiseaseParameters",
    "EpidemicModel",
    "InterventionPlan",
    "intervene_bridge_links",
    "intervene_within_groups",
    "PopulationParameters",
    "preset_general_cluster",
    "RK4Stepper",
    "SimulationParameters",
    "StatePoint",
    "StateStepper",
    "StructuredEpidemicModel",
    "StructuredStatePoint",
]
__version__ = "0.1.0"
