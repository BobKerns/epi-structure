"""Parameter definitions for epidemic population structure models."""

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class SimulationParameters:
    """Global simulation controls."""

    time_step: float = 0.1
    duration: float = 100.0


@dataclass(frozen=True, slots=True)
class PopulationParameters:
    """Per-population epidemiological parameters."""

    size: float
    beta: float
    latent_period: float | None = None
    infectious_period: float = 7.0
    waning_period: float | None = None
