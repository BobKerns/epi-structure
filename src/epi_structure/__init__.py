"""Epidemic population structure simulation package."""

from .contact_matrix import build_contact_matrix, preset_general_cluster
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
    "DiseaseParameters",
    "EpidemicModel",
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
