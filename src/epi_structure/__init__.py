"""Epidemic population structure simulation package."""

from .model import EpidemicModel, RK4Stepper, StatePoint, StateStepper
from .parameters import DiseaseParameters, PopulationParameters, SimulationParameters

__all__ = [
    "DiseaseParameters",
    "EpidemicModel",
    "PopulationParameters",
    "RK4Stepper",
    "SimulationParameters",
    "StatePoint",
    "StateStepper",
]
__version__ = "0.1.0"
