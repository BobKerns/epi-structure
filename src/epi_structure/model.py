"""Core epidemic model implementation."""

import math
import random
from dataclasses import dataclass
from typing import Callable, Protocol

import pandas as pd

from .interventions import InterventionPlan
from .parameters import PopulationParameters, SimulationParameters


StateVector = tuple[float, ...]
DerivativeFunction = Callable[[StateVector], StateVector]


class StateStepper(Protocol):
    """Strategy interface for advancing state by one simulation step."""

    def step(
        self,
        state: StateVector,
        time_step: float,
        derivatives: DerivativeFunction,
        context: object | None = None,
    ) -> StateVector:
        """Advance one integration/transition step and return the next state."""


class RK4Stepper:
    """Classical fourth-order Runge-Kutta state stepper."""

    @staticmethod
    def _add_scaled(base: StateVector, delta: StateVector, scale: float) -> StateVector:
        return tuple(base[i] + delta[i] * scale for i in range(len(base)))

    def step(
        self,
        state: StateVector,
        time_step: float,
        derivatives: DerivativeFunction,
        context: object | None = None,
    ) -> StateVector:
        k1 = derivatives(state)
        k2 = derivatives(self._add_scaled(state, k1, 0.5 * time_step))
        k3 = derivatives(self._add_scaled(state, k2, 0.5 * time_step))
        k4 = derivatives(self._add_scaled(state, k3, time_step))

        return tuple(
            state[i] + time_step * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
            for i in range(len(state))
        )


class TransitionProbabilityStepper:
    """Discrete-time stochastic transition engine.

    Each compartment transition is sampled as a binomial draw from the current
    compartment counts, so the state remains whole-numbered rather than
    containing fractional people.
    """

    def __init__(self, seed: int | None = None) -> None:
        self._rng = random.Random(seed)

    @staticmethod
    def _transition_probability(rate: float, time_step: float) -> float:
        if rate <= 0.0 or time_step <= 0.0:
            return 0.0
        return min(1.0, -math.expm1(-rate * time_step))

    def _binomial(self, trials: float, probability: float) -> int:
        count = max(0, int(round(trials)))
        probability = min(max(probability, 0.0), 1.0)
        if count == 0 or probability == 0.0:
            return 0
        if probability == 1.0:
            return count
        return sum(1 for _ in range(count) if self._rng.random() < probability)

    def _step_single_population(self, context, state: StateVector, time_step: float) -> StateVector:
        population = context.population
        compartments = population.disease.compartments
        has_deceased = "D" in compartments

        if has_deceased:
            if len(state) != 5:
                raise ValueError(f"SEIRD model requires 5-element state, got {len(state)}")
            susceptible, exposed, infected, recovered, deceased = state
        else:
            if len(state) != 4:
                raise ValueError(f"SEIR model requires 4-element state, got {len(state)}")
            susceptible, exposed, infected, recovered = state
            deceased = 0.0

        susceptible = max(0.0, susceptible)
        exposed = max(0.0, exposed)
        infected = max(0.0, infected)
        recovered = max(0.0, recovered)
        deceased = max(0.0, deceased)

        infection_rate = 0.0 if population.size <= 0 else population.beta * infected / float(population.size)
        infection_probability = self._transition_probability(infection_rate, time_step)
        waning_probability = self._transition_probability(context._omega, time_step)

        if context._sigma is None:
            new_infections = self._binomial(susceptible, infection_probability)
            exposed_to_infected = 0
        else:
            new_infections = self._binomial(susceptible, infection_probability)
            exposed_to_infected = self._binomial(
                exposed,
                self._transition_probability(context._sigma, time_step),
            )

        infected_leavers = self._binomial(
            infected,
            self._transition_probability(context._gamma, time_step),
        )

        if has_deceased:
            cfr = population.disease.case_fatality_rate
            deaths = self._binomial(infected_leavers, cfr)
            recoveries = infected_leavers - deaths
        else:
            deaths = 0
            recoveries = infected_leavers

        waning = self._binomial(recovered, waning_probability)

        if context._sigma is None:
            next_susceptible = susceptible - new_infections + waning
            next_exposed = 0.0
            next_infected = infected + new_infections - infected_leavers
            next_recovered = recovered + recoveries - waning
        else:
            next_susceptible = susceptible - new_infections + waning
            next_exposed = exposed + new_infections - exposed_to_infected
            next_infected = infected + exposed_to_infected - infected_leavers
            next_recovered = recovered + recoveries - waning

        if has_deceased:
            next_deceased = deceased + deaths
            return (
                float(next_susceptible),
                float(next_exposed),
                float(next_infected),
                float(next_recovered),
                float(next_deceased),
            )

        return (
            float(next_susceptible),
            float(next_exposed),
            float(next_infected),
            float(next_recovered),
        )

    def _step_structured(self, context, state: StateVector, time_step: float) -> StateVector:
        compartments = context._compartments
        offsets = context._compartment_offsets
        has_deceased = context._has_deceased
        next_state = list(state)

        infected = [
            max(0.0, state[context._state_index(i, offsets["I"])])
            for i in range(len(context.populations))
        ]

        for i in range(len(context.populations)):
            s_index = context._state_index(i, offsets["S"])
            i_index = context._state_index(i, offsets["I"])
            r_index = context._state_index(i, offsets["R"])
            e_index = None if "E" not in offsets else context._state_index(i, offsets["E"])
            d_index = None if "D" not in offsets else context._state_index(i, offsets["D"])

            susceptible = max(0.0, state[s_index])
            exposed = 0.0 if e_index is None else max(0.0, state[e_index])
            infected_count = max(0.0, state[i_index])
            recovered = max(0.0, state[r_index])
            deceased = 0.0 if d_index is None else max(0.0, state[d_index])

            infection_pressure = sum(
                context._active_contact_matrix[i][j] * (infected[j] / context._sizes[j])
                for j in range(len(context.populations))
            )
            infection_probability = self._transition_probability(infection_pressure, time_step)
            waning_probability = self._transition_probability(context._omegas[i], time_step)

            new_infections = self._binomial(susceptible, infection_probability)
            if e_index is None:
                exposed_to_infected = 0
            else:
                exposed_to_infected = self._binomial(
                    exposed,
                    self._transition_probability(context._sigmas[i], time_step),
                )

            infected_leavers = self._binomial(
                infected_count,
                self._transition_probability(context._gammas[i], time_step),
            )

            if has_deceased:
                deaths = self._binomial(infected_leavers, context._cfrs[i])
                recoveries = infected_leavers - deaths
            else:
                deaths = 0
                recoveries = infected_leavers

            waning = self._binomial(recovered, waning_probability)

            next_state[s_index] = float(susceptible - new_infections + waning)
            if e_index is not None:
                next_state[e_index] = float(exposed + new_infections - exposed_to_infected)
            next_state[i_index] = float(infected_count + exposed_to_infected - infected_leavers)
            next_state[r_index] = float(recovered + recoveries - waning)
            if d_index is not None:
                next_state[d_index] = float(deceased + deaths)

        return tuple(next_state)

    def step(
        self,
        state: StateVector,
        time_step: float,
        derivatives: DerivativeFunction,
        context: object | None = None,
    ) -> StateVector:
        if context is None:
            raise ValueError("TransitionProbabilityStepper requires model context")

        if hasattr(context, "populations") and hasattr(context, "_compartment_offsets"):
            return self._step_structured(context, state, time_step)

        return self._step_single_population(context, state, time_step)


def _step_with_context(
    stepper: StateStepper,
    state: StateVector,
    time_step: float,
    derivatives: DerivativeFunction,
    context: object,
) -> StateVector:
    try:
        return stepper.step(
            state=state,
            time_step=time_step,
            derivatives=derivatives,
            context=context,
        )
    except TypeError as error:
        if "context" not in str(error):
            raise
        return stepper.step(state=state, time_step=time_step, derivatives=derivatives)


@dataclass(frozen=True, slots=True)
class StatePoint:
    """State snapshot at a single simulation time.
    
    Legacy SEIR-specific format. Maintained for backward compatibility.
    For generalized compartments, use GeneralizedStatePoint.
    """

    time: float
    susceptible: float
    exposed: float
    infected: float
    recovered: float

    @property
    def total_population(self) -> float:
        """Return total population represented by this state."""

        return self.susceptible + self.exposed + self.infected + self.recovered


@dataclass(frozen=True, slots=True)
class GeneralizedStatePoint:
    """State snapshot for flexible compartment models (SEIR, SEIRD, etc.).
    
    Stores compartment values in a dict, keyed by compartment name.
    Supports arbitrary compartment lists while maintaining a clear interface.
    """
    
    time: float
    compartments: tuple[str, ...]  # Ordered list of compartment names
    values: dict[str, float]  # Compartment name -> current count/proportion
    
    @property
    def total_population(self) -> float:
        """Return sum of all compartments (population size)."""
        return sum(self.values.values())
    
    def get(self, compartment: str) -> float:
        """Get value for a specific compartment."""
        return self.values.get(compartment, 0.0)


@dataclass(frozen=True, slots=True)
class StructuredStatePoint:
    """State snapshot for a multi-population simulation timepoint."""

    time: float
    by_population: dict[str, StatePoint]

    @property
    def total_population(self) -> float:
        """Return total population represented across all populations."""

        return sum(point.total_population for point in self.by_population.values())


@dataclass(slots=True)
class ScenarioState:
    """Runtime simulation state for a multi-population scenario.

    ScenarioState holds the current state of a scenario at a point in time.
    It is initialized from a Scenario (immutable config) and evolves through
    simulation steps. The Model class handles state transitions; ScenarioState
    manages persistence and by-population snapshots.

    Attributes:
        time: Current simulation time.
        by_population: Mapping of population name -> current StatePoint.
        active_contact_matrix: Current contact matrix (may change via interventions).
        population_names: List of population names (from scenario config).
        population_sizes: List of population sizes (for normalization).
    """

    time: float
    by_population: dict[str, StatePoint]
    active_contact_matrix: list[list[float]]
    population_names: list[str]
    population_sizes: list[float]

    @classmethod
    def from_scenario(cls, scenario) -> "ScenarioState":  # noqa: F821
        """Initialize ScenarioState from a Scenario configuration.

        Args:
            scenario: A Scenario object with populations, contact_matrix, etc.

        Returns:
            ScenarioState at t=0 with initial populations from scenario.
        """
        from .scenarios import Scenario  # Import here to avoid circular dependency

        by_population = {
            pop.name: StatePoint(
                time=0.0,
                susceptible=float(pop.initial_susceptible),
                exposed=float(pop.initial_exposed),
                infected=float(pop.initial_infected),
                recovered=float(pop.initial_recovered),
            )
            for pop in scenario.populations
        }

        population_sizes = [float(pop.size) for pop in scenario.populations]

        return cls(
            time=0.0,
            by_population=by_population,
            active_contact_matrix=[row[:] for row in scenario.contact_matrix],
            population_names=list(scenario.population_names),
            population_sizes=population_sizes,
        )

    @property
    def total_population(self) -> float:
        """Return total population across all groups."""
        return sum(state.total_population for state in self.by_population.values())

    @property
    def state_vector(self) -> StateVector:
        """Extract current state as a flat tuple for numerical integration.

        Order: [S_0, E_0, I_0, R_0, S_1, E_1, I_1, R_1, ...]
        """
        values: list[float] = []
        for name in self.population_names:
            state = self.by_population[name]
            values.extend((state.susceptible, state.exposed, state.infected, state.recovered))
        return tuple(values)

    def update_from_state_vector(self, state_vector: StateVector, time: float) -> None:
        """Update by_population and time from a state vector.

        Used after numerical integration step to update the scenario state.
        Also normalizes compartments to ensure exact population conservation.

        Args:
            state_vector: Flat tuple [S_0, E_0, I_0, R_0, ...] from integration.
            time: New simulation time.
        """
        self.time = time
        offset = 0
        for i, name in enumerate(self.population_names):
            s = max(0.0, state_vector[offset])
            e = max(0.0, state_vector[offset + 1])
            infected = max(0.0, state_vector[offset + 2])
            recovered = max(0.0, state_vector[offset + 3])

            total = s + e + infected + recovered
            correction = self.population_sizes[i] - total
            s = max(0.0, s + correction)

            self.by_population[name] = StatePoint(
                time=time,
                susceptible=s,
                exposed=e,
                infected=infected,
                recovered=recovered,
            )
            offset += 4


class EpidemicModel:
    """Deterministic single-population SEIRS simulator.

    The time-stepping mechanism is pluggable through the ``stepper`` argument.
    ``RK4Stepper`` is the default for current Phase 1 behavior.
    """

    def __init__(
        self,
        population: PopulationParameters,
        simulation: SimulationParameters,
        stepper: StateStepper | None = None,
    ) -> None:
        self.population = population
        self.simulation = simulation
        self.stepper = stepper or RK4Stepper()

        self._gamma = 1.0 / self.population.disease.infectious_period
        self._sigma = (
            None
            if self.population.disease.latent_period is None
            else 1.0 / self.population.disease.latent_period
        )
        self._omega = (
            0.0
            if self.population.disease.waning_period is None
            else 1.0 / self.population.disease.waning_period
        )

    def _derivatives(self, state: StateVector) -> StateVector:
        """Compute derivatives for current state.
        
        Supports both SEIR (4 compartments) and SEIRD (5 compartments).
        For SEIRD, the case_fatality_rate determines the split between
        recovery (1-CFR) and death (CFR) flows from infectious state.
        """
        compartments = self.population.disease.compartments
        has_deceased = "D" in compartments
        
        if has_deceased:
            # SEIRD model: S, E, I, R, D
            if len(state) != 5:
                raise ValueError(f"SEIRD model requires 5-element state, got {len(state)}")
            s, e, i, r, d = state
        else:
            # SEIR model: S, E, I, R
            if len(state) != 4:
                raise ValueError(f"SEIR model requires 4-element state, got {len(state)}")
            s, e, i, r = state
            d = None
        
        infection_flow = self.population.beta * s * i / self.population.size
        
        # Base transitions (S → E → I → R/D with optional waning)
        if self._sigma is None:
            # SIR variant: no exposed state
            ds_dt = -infection_flow + self._omega * r
            de_dt = 0.0
            di_dt = infection_flow - self._gamma * i
            dr_dt = self._gamma * i - self._omega * r
            if has_deceased:
                # Split recovery into R and D based on CFR
                cfr = self.population.disease.case_fatality_rate
                dr_dt = (1.0 - cfr) * self._gamma * i - self._omega * r
                dd_dt = cfr * self._gamma * i
        else:
            # SEIR(D) variant: exposed state present
            ds_dt = -infection_flow + self._omega * r
            de_dt = infection_flow - self._sigma * e
            di_dt = self._sigma * e - self._gamma * i
            dr_dt = self._gamma * i - self._omega * r
            if has_deceased:
                # Split recovery into R and D based on CFR
                cfr = self.population.disease.case_fatality_rate
                dr_dt = (1.0 - cfr) * self._gamma * i - self._omega * r
                dd_dt = cfr * self._gamma * i
        
        if has_deceased:
            return ds_dt, de_dt, di_dt, dr_dt, dd_dt
        else:
            return ds_dt, de_dt, di_dt, dr_dt

    @staticmethod
    def _normalize_state(
        state: StateVector,
        expected_total: float,
        compartments: list[str],
    ) -> StateVector:
        """Normalize state vector to ensure exact population conservation.
        
        Guards against small negative drift from finite-precision arithmetic.
        Works for both SEIR (4 elements) and SEIRD (5 elements).
        """
        has_deceased = "D" in compartments
        
        # Guard against small negative drift
        if has_deceased:
            s, e, i, r, d = state
            s = max(0.0, s)
            e = max(0.0, e)
            i = max(0.0, i)
            r = max(0.0, r)
            d = max(0.0, d)
            total = s + e + i + r + d
            correction = expected_total - total
            s = max(0.0, s + correction)
            return s, e, i, r, d
        else:
            s = max(0.0, state[0])
            e = max(0.0, state[1])
            i = max(0.0, state[2])
            r = max(0.0, state[3])
            total = s + e + i + r
            correction = expected_total - total
            s = max(0.0, s + correction)
            return s, e, i, r

    def simulate(self, tidy: bool = False):
        """Run the simulation and return sampled state points or a tidy DataFrame if tidy=True.
        
        Supports both SEIR and SEIRD models based on disease.compartments.
        """

        dt = self.simulation.time_step
        steps = int(round(self.simulation.duration / dt))
        
        # Build initial state vector from compartments in order
        compartments = self.population.disease.compartments
        has_deceased = "D" in compartments
        initial_state_dict = self.population.initial_state_by_compartment()
        state: StateVector = tuple(initial_state_dict[c] for c in compartments)

        # Keep raw state snapshots so tidy output can include all compartments.
        trajectory: list[tuple[float, StateVector]] = [(0.0, state)]

        for step in range(1, steps + 1):
            raw_state = _step_with_context(self.stepper, state, dt, self._derivatives, self)
            state = self._normalize_state(raw_state, float(self.population.size), compartments)
            if step % self.simulation.output_stride == 0:
                trajectory.append((step * dt, state))

        if not tidy:
            if not has_deceased:
                return [StatePoint(t, s[0], s[1], s[2], s[3]) for t, s in trajectory]

            generalized: list[GeneralizedStatePoint] = []
            for t, s in trajectory:
                values = {compartment: s[i] for i, compartment in enumerate(compartments)}
                generalized.append(
                    GeneralizedStatePoint(
                        time=t,
                        compartments=tuple(compartments),
                        values=values,
                    )
                )
            return generalized

        # Tidy DataFrame: one row per time, one column per compartment
        import pandas as pd

        records: list[dict[str, float]] = []
        for t, s in trajectory:
            row: dict[str, float] = {
                "time": t,
                "susceptible": s[compartments.index("S")] if "S" in compartments else 0.0,
                "exposed": s[compartments.index("E")] if "E" in compartments else 0.0,
                "infected": s[compartments.index("I")] if "I" in compartments else 0.0,
                "recovered": s[compartments.index("R")] if "R" in compartments else 0.0,
                "total_population": sum(s),
            }
            if has_deceased:
                row["deceased"] = s[compartments.index("D")]
            records.append(row)
        return pd.DataFrame(records)


class StructuredEpidemicModel:
    """Deterministic multi-population SEIRS simulator with contact matrix coupling.

    The force of infection on population i is:
        lambda_i = S_i * sum_j(beta_ij * I_j)
    where beta_ij is the contact matrix entry at row i, column j.
    """

    def __init__(
        self,
        populations: list[PopulationParameters],
        contact_matrix: list[list[float]],
        simulation: SimulationParameters,
        intervention_plan: InterventionPlan | None = None,
        stepper: StateStepper | None = None,
    ) -> None:
        if not populations:
            raise ValueError("populations must not be empty")

        self.populations = populations
        self.simulation = simulation
        self.stepper = stepper or RK4Stepper()
        self.intervention_plan = intervention_plan
        self.contact_matrix = self._validate_contact_matrix(contact_matrix, len(populations))

        self._names = [population.name for population in populations]
        self._compartments = list(populations[0].disease.compartments)
        for population in populations:
            if list(population.disease.compartments) != self._compartments:
                raise ValueError("all populations must use the same compartment ordering")
        self._compartment_count = len(self._compartments)
        self._compartment_offsets = {
            compartment: index for index, compartment in enumerate(self._compartments)
        }
        self._has_deceased = "D" in self._compartments

        required = {"S", "I", "R"}
        if not required.issubset(self._compartment_offsets):
            raise ValueError("compartments must include S, I, and R for StructuredEpidemicModel")

        self._sizes = [float(population.size) for population in populations]
        self._gammas = [1.0 / population.disease.infectious_period for population in populations]
        self._sigmas = [
            None if population.disease.latent_period is None else 1.0 / population.disease.latent_period
            for population in populations
        ]
        if any(sigma is not None for sigma in self._sigmas) and "E" not in self._compartment_offsets:
            raise ValueError("compartments must include E when latent_period is provided")
        self._omegas = [
            0.0 if population.disease.waning_period is None else 1.0 / population.disease.waning_period
            for population in populations
        ]
        self._cfrs = [population.disease.case_fatality_rate for population in populations]
        self._active_contact_matrix = [row[:] for row in self.contact_matrix]

    @classmethod
    def from_scenario(cls, scenario) -> "StructuredEpidemicModel":  # noqa: F821
        """Create a StructuredEpidemicModel from a Scenario configuration.

        Args:
            scenario: A Scenario object with populations, contact_matrix, etc.

        Returns:
            StructuredEpidemicModel configured from the scenario.
        """
        return cls(
            populations=scenario.populations,
            contact_matrix=scenario.contact_matrix,
            simulation=scenario.simulation,
            intervention_plan=scenario.intervention_plan,
        )

    @staticmethod
    def _validate_contact_matrix(contact_matrix: list[list[float]], count: int) -> list[list[float]]:
        if len(contact_matrix) != count:
            raise ValueError("contact_matrix must have one row per population")
        for row in contact_matrix:
            if len(row) != count:
                raise ValueError("contact_matrix must be square with one column per population")
            if any(value < 0 for value in row):
                raise ValueError("contact_matrix entries must be non-negative")
        return contact_matrix

    def _state_index(self, population_index: int, compartment_offset: int) -> int:
        return self._compartment_count * population_index + compartment_offset

    def _initial_state(self) -> StateVector:
        values: list[float] = []
        for population in self.populations:
            initial_state = population.initial_state_by_compartment()
            for compartment in self._compartments:
                values.append(initial_state.get(compartment, 0.0))
        return tuple(values)

    def _derivatives(self, state: StateVector) -> StateVector:
        derivatives: list[float] = [0.0] * len(state)
        infected_offset = self._compartment_offsets["I"]
        infected = [state[self._state_index(i, infected_offset)] for i in range(len(self.populations))]

        for i in range(len(self.populations)):
            s_index = self._state_index(i, self._compartment_offsets["S"])
            i_index = self._state_index(i, self._compartment_offsets["I"])
            r_index = self._state_index(i, self._compartment_offsets["R"])
            e_index = (
                None
                if "E" not in self._compartment_offsets
                else self._state_index(i, self._compartment_offsets["E"])
            )
            d_index = (
                None
                if "D" not in self._compartment_offsets
                else self._state_index(i, self._compartment_offsets["D"])
            )

            s_value = state[s_index]
            e_value = 0.0 if e_index is None else state[e_index]
            i_value = state[i_index]
            r_value = state[r_index]

            infection_pressure = sum(
                self._active_contact_matrix[i][j] * (infected[j] / self._sizes[j])
                for j in range(len(self.populations))
            )
            infection_flow = s_value * infection_pressure

            sigma = self._sigmas[i]
            gamma = self._gammas[i]
            omega = self._omegas[i]

            if sigma is None:
                ds_dt = -infection_flow + omega * r_value
                de_dt = 0.0
                di_dt = infection_flow - gamma * i_value
                dr_dt = gamma * i_value - omega * r_value
            else:
                ds_dt = -infection_flow + omega * r_value
                de_dt = infection_flow - sigma * e_value
                di_dt = sigma * e_value - gamma * i_value
                dr_dt = gamma * i_value - omega * r_value

            if d_index is not None:
                cfr = self._cfrs[i]
                dr_dt = (1.0 - cfr) * gamma * i_value - omega * r_value
                derivatives[d_index] = cfr * gamma * i_value

            derivatives[s_index] = ds_dt
            if e_index is not None:
                derivatives[e_index] = de_dt
            derivatives[i_index] = di_dt
            derivatives[r_index] = dr_dt

        return tuple(derivatives)

    def _normalize_state(self, state: StateVector) -> StateVector:
        normalized: list[float] = list(state)

        for i in range(len(self.populations)):
            s_index = self._state_index(i, self._compartment_offsets["S"])

            compartment_indices = [
                self._state_index(i, offset)
                for offset in self._compartment_offsets.values()
            ]
            for index in compartment_indices:
                normalized[index] = max(0.0, normalized[index])

            total = sum(normalized[index] for index in compartment_indices)
            correction = self._sizes[i] - total
            normalized[s_index] = max(0.0, normalized[s_index] + correction)

        return tuple(normalized)

    def _structured_state_point(self, time: float, state: StateVector) -> StructuredStatePoint:
        by_population: dict[str, StatePoint] = {}
        for i, name in enumerate(self._names):
            offset = self._state_index(i, 0)
            s = state[offset + self._compartment_offsets["S"]]
            e = (
                0.0
                if "E" not in self._compartment_offsets
                else state[offset + self._compartment_offsets["E"]]
            )
            infected = state[offset + self._compartment_offsets["I"]]
            recovered = state[offset + self._compartment_offsets["R"]]
            by_population[name] = StatePoint(
                time=time,
                susceptible=s,
                exposed=e,
                infected=infected,
                recovered=recovered,
            )
        return StructuredStatePoint(time=time, by_population=by_population)

    def simulate(self, tidy: bool = False):
        """Run the simulation and return sampled multi-population state points or a tidy DataFrame if tidy=True.

        The simulation uses ScenarioState internally to track runtime state evolution,
        providing a clear separation between configuration (scenario) and runtime state.

        Args:
            tidy: If False (default), return list of StructuredStatePoint.
                 If True, return tidy pandas DataFrame with one row per (time, population).

        Returns:
            List[StructuredStatePoint] or pd.DataFrame depending on tidy parameter.
        """

        dt = self.simulation.time_step
        steps = int(round(self.simulation.duration / dt))

        # Initialize scenario state from populations (mimics ScenarioState initialization)
        current_state = self._initial_state()

        trajectory = [self._structured_state_point(0.0, current_state)]

        for step in range(1, steps + 1):
            current_time = (step - 1) * dt
            if self.intervention_plan is None:
                self._active_contact_matrix = [row[:] for row in self.contact_matrix]
            else:
                self._active_contact_matrix = self.intervention_plan.matrix_at_time(
                    base_contact_matrix=self.contact_matrix,
                    population_names=self._names,
                    time=current_time,
                )

            raw_state = _step_with_context(
                self.stepper,
                current_state,
                dt,
                self._derivatives,
                self,
            )
            current_state = self._normalize_state(raw_state)
            if step % self.simulation.output_stride == 0:
                trajectory.append(self._structured_state_point(step * dt, current_state))

        if not tidy:
            return trajectory

        # Tidy DataFrame: one row per time, population, compartment
        import pandas as pd
        records = []
        for pt in trajectory:
            for pop_name, pop_state in pt.by_population.items():
                record = {
                    "time": pt.time,
                    "population": pop_name,
                    "susceptible": pop_state.susceptible,
                    "exposed": pop_state.exposed,
                    "infected": pop_state.infected,
                    "recovered": pop_state.recovered,
                    "total_population": pop_state.total_population,
                }
                if self._has_deceased:
                    # Deceased is represented in runtime state but not in StatePoint,
                    # so infer it from the population size and visible compartments.
                    size = next(p.size for p in self.populations if p.name == pop_name)
                    record["deceased"] = max(0.0, size - (
                        pop_state.susceptible + pop_state.exposed + pop_state.infected + pop_state.recovered
                    ))
                records.append(record)
        return pd.DataFrame(records)

    def contact_matrix_df(self) -> pd.DataFrame:
        """Return the contact matrix as a pandas DataFrame with population names as index/columns."""
        return pd.DataFrame(self.contact_matrix, index=self._names, columns=self._names)
