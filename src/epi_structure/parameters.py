"""Parameter definitions for epidemic population structure models."""

from dataclasses import dataclass, field


@dataclass(frozen=True, slots=True)
class SimulationParameters:
    """Global simulation controls.

    Attributes:
        time_step: Integration step size in time units.
        duration: Total simulated time span in time units.
        output_stride: Emit every Nth integration step.
    """

    time_step: float = 0.1
    duration: float = 100.0
    output_stride: int = 1

    def __post_init__(self) -> None:
        if self.time_step <= 0:
            raise ValueError("time_step must be positive")
        if self.duration <= 0:
            raise ValueError("duration must be positive")
        if self.output_stride < 1:
            raise ValueError("output_stride must be at least 1")


@dataclass(frozen=True, slots=True)
class DiseaseParameters:
    """Shared progression parameters for a single population or group.
    
    Supports generalized compartment models (SEIR, SEIRD, SEIARD, etc.) via
    the compartments field. Mortality is modeled via case_fatality_rate (CFR),
    which diverts individuals from recovery to a death compartment.
    """

    infectious_period: float = 7.0
    latent_period: float | None = None
    waning_period: float | None = None
    compartments: list[str] = field(default_factory=lambda: ["S", "E", "I", "R"])
    case_fatality_rate: float = 0.0

    def __post_init__(self) -> None:
        if self.infectious_period <= 0:
            raise ValueError("infectious_period must be positive")
        if self.latent_period is not None and self.latent_period <= 0:
            raise ValueError("latent_period must be positive when provided")
        if self.waning_period is not None and self.waning_period <= 0:
            raise ValueError("waning_period must be positive when provided")
        
        # Validate compartments list
        if not self.compartments:
            raise ValueError("compartments list must not be empty")
        if len(set(self.compartments)) != len(self.compartments):
            raise ValueError("compartments list has duplicates")
        if "S" not in self.compartments:
            raise ValueError("compartments must include 'S' (Susceptible)")
        
        # Validate CFR
        if not 0.0 <= self.case_fatality_rate <= 1.0:
            raise ValueError("case_fatality_rate must be in [0.0, 1.0]")
        
        # If D (deceased) compartment exists, CFR must be > 0
        if "D" in self.compartments and self.case_fatality_rate == 0.0:
            raise ValueError("compartments include 'D' but case_fatality_rate is 0")
        
        # If CFR > 0, D compartment must be present
        if self.case_fatality_rate > 0.0 and "D" not in self.compartments:
            raise ValueError("case_fatality_rate > 0 but 'D' (Deceased) not in compartments")


@dataclass(frozen=True, slots=True)
class PopulationParameters:
    """Per-population state and epidemiological parameters.

    Attributes:
        name: Human-readable population label.
        size: Total number of individuals in the population.
        beta: Within-population transmission coefficient.
        initial_susceptible: Susceptible individuals at t=0.
        initial_exposed: Exposed individuals at t=0.
        initial_infected: Infectious individuals at t=0.
        initial_recovered: Recovered individuals at t=0.
        initial_deceased: Deceased individuals at t=0 (for SEIRD models).
        disease: Progression parameters for this population.
    """

    name: str
    size: int
    beta: float
    initial_susceptible: int | None = None
    initial_exposed: int = 0
    initial_infected: int = 1
    initial_recovered: int = 0
    initial_deceased: int = 0
    disease: DiseaseParameters = field(default_factory=DiseaseParameters)

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("name must not be empty")
        if self.size <= 0:
            raise ValueError("size must be positive")
        if self.beta < 0:
            raise ValueError("beta must be non-negative")
        if self.initial_exposed < 0:
            raise ValueError("initial_exposed must be non-negative")
        if self.initial_infected < 0:
            raise ValueError("initial_infected must be non-negative")
        if self.initial_recovered < 0:
            raise ValueError("initial_recovered must be non-negative")
        if self.initial_deceased < 0:
            raise ValueError("initial_deceased must be non-negative")

        if self.initial_susceptible is None:
            susceptible = self.size - self.initial_exposed - self.initial_infected - self.initial_recovered - self.initial_deceased
            if susceptible < 0:
                raise ValueError("initial compartment counts exceed population size")
            object.__setattr__(self, "initial_susceptible", susceptible)
        elif self.initial_susceptible < 0:
            raise ValueError("initial_susceptible must be non-negative")

        total_initial = (
            self.initial_susceptible
            + self.initial_exposed
            + self.initial_infected
            + self.initial_recovered
            + self.initial_deceased
        )
        if total_initial != self.size:
            raise ValueError("initial compartment counts must sum to population size")

    def initial_state_by_compartment(self) -> dict[str, float]:
        """Build initial state dict keyed by compartment name.
        
        Maps compartment names from disease.compartments to initial values.
        Supports both SEIR (hardcoded fields) and SEIRD+ (with D) seamlessly.
        """
        state_dict = {}
        for compartment in self.disease.compartments:
            if compartment == "S":
                state_dict[compartment] = float(self.initial_susceptible)
            elif compartment == "E":
                state_dict[compartment] = float(self.initial_exposed)
            elif compartment == "I":
                state_dict[compartment] = float(self.initial_infected)
            elif compartment == "R":
                state_dict[compartment] = float(self.initial_recovered)
            elif compartment == "D":
                state_dict[compartment] = float(self.initial_deceased)
            else:
                # Unknown compartment; default to 0
                state_dict[compartment] = 0.0
        return state_dict

    @property
    def initial_state(self) -> tuple[int, int, int, int]:
        """Return the initial S, E, I, R state as a tuple."""

        return (
            self.initial_susceptible,
            self.initial_exposed,
            self.initial_infected,
            self.initial_recovered,
        )
