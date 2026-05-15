"""Core epidemic model implementation."""

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
    ) -> StateVector:
        k1 = derivatives(state)
        k2 = derivatives(self._add_scaled(state, k1, 0.5 * time_step))
        k3 = derivatives(self._add_scaled(state, k2, 0.5 * time_step))
        k4 = derivatives(self._add_scaled(state, k3, time_step))

        return tuple(
            state[i] + time_step * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
            for i in range(len(state))
        )


@dataclass(frozen=True, slots=True)
class StatePoint:
    """State snapshot at a single simulation time."""

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
class StructuredStatePoint:
    """State snapshot for a multi-population simulation timepoint."""

    time: float
    by_population: dict[str, StatePoint]

    @property
    def total_population(self) -> float:
        """Return total population represented across all populations."""

        return sum(point.total_population for point in self.by_population.values())


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
        s, e, i, r = state
        infection_flow = self.population.beta * s * i

        if self._sigma is None:
            ds_dt = -infection_flow + self._omega * r
            de_dt = 0.0
            di_dt = infection_flow - self._gamma * i
            dr_dt = self._gamma * i - self._omega * r
        else:
            ds_dt = -infection_flow + self._omega * r
            de_dt = infection_flow - self._sigma * e
            di_dt = self._sigma * e - self._gamma * i
            dr_dt = self._gamma * i - self._omega * r

        return ds_dt, de_dt, di_dt, dr_dt

    @staticmethod
    def _normalize_state(state: StateVector, expected_total: float) -> StateVector:
        # Guard against small negative drift from finite-precision arithmetic.
        s = max(0.0, state[0])
        e = max(0.0, state[1])
        i = max(0.0, state[2])
        r = max(0.0, state[3])
        total = s + e + i + r
        correction = expected_total - total
        s = max(0.0, s + correction)
        return s, e, i, r

    def simulate(self, tidy: bool = False):
        """Run the simulation and return sampled state points or a tidy DataFrame if tidy=True."""

        dt = self.simulation.time_step
        steps = int(round(self.simulation.duration / dt))
        state: StateVector = tuple(float(v) for v in self.population.initial_state)

        trajectory = [StatePoint(0.0, *state)]

        for step in range(1, steps + 1):
            raw_state = self.stepper.step(state=state, time_step=dt, derivatives=self._derivatives)
            state = self._normalize_state(raw_state, float(self.population.size))
            if step % self.simulation.output_stride == 0:
                trajectory.append(StatePoint(step * dt, *state))

        if not tidy:
            return trajectory

        # Tidy DataFrame: one row per time, one column per compartment
        import pandas as pd
        records = [
            {
                "time": pt.time,
                "susceptible": pt.susceptible,
                "exposed": pt.exposed,
                "infected": pt.infected,
                "recovered": pt.recovered,
                "total_population": pt.total_population,
            }
            for pt in trajectory
        ]
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
        self._sizes = [float(population.size) for population in populations]
        self._gammas = [1.0 / population.disease.infectious_period for population in populations]
        self._sigmas = [
            None if population.disease.latent_period is None else 1.0 / population.disease.latent_period
            for population in populations
        ]
        self._omegas = [
            0.0 if population.disease.waning_period is None else 1.0 / population.disease.waning_period
            for population in populations
        ]
        self._active_contact_matrix = [row[:] for row in self.contact_matrix]

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

    @staticmethod
    def _state_index(population_index: int, compartment_offset: int) -> int:
        return 4 * population_index + compartment_offset

    def _initial_state(self) -> StateVector:
        values: list[float] = []
        for population in self.populations:
            s, e, i, r = population.initial_state
            values.extend((float(s), float(e), float(i), float(r)))
        return tuple(values)

    def _derivatives(self, state: StateVector) -> StateVector:
        derivatives: list[float] = [0.0] * len(state)
        infected = [state[self._state_index(i, 2)] for i in range(len(self.populations))]

        for i in range(len(self.populations)):
            s_index = self._state_index(i, 0)
            e_index = self._state_index(i, 1)
            i_index = self._state_index(i, 2)
            r_index = self._state_index(i, 3)

            s_value = state[s_index]
            e_value = state[e_index]
            i_value = state[i_index]
            r_value = state[r_index]

            infection_pressure = sum(
                self._active_contact_matrix[i][j] * infected[j] for j in range(len(self.populations))
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

            derivatives[s_index] = ds_dt
            derivatives[e_index] = de_dt
            derivatives[i_index] = di_dt
            derivatives[r_index] = dr_dt

        return tuple(derivatives)

    def _normalize_state(self, state: StateVector) -> StateVector:
        normalized: list[float] = list(state)

        for i in range(len(self.populations)):
            s_index = self._state_index(i, 0)
            e_index = self._state_index(i, 1)
            i_index = self._state_index(i, 2)
            r_index = self._state_index(i, 3)

            s = max(0.0, normalized[s_index])
            e = max(0.0, normalized[e_index])
            infected = max(0.0, normalized[i_index])
            recovered = max(0.0, normalized[r_index])

            total = s + e + infected + recovered
            correction = self._sizes[i] - total
            s = max(0.0, s + correction)

            normalized[s_index] = s
            normalized[e_index] = e
            normalized[i_index] = infected
            normalized[r_index] = recovered

        return tuple(normalized)

    def _structured_state_point(self, time: float, state: StateVector) -> StructuredStatePoint:
        by_population: dict[str, StatePoint] = {}
        for i, name in enumerate(self._names):
            offset = self._state_index(i, 0)
            by_population[name] = StatePoint(
                time=time,
                susceptible=state[offset],
                exposed=state[offset + 1],
                infected=state[offset + 2],
                recovered=state[offset + 3],
            )
        return StructuredStatePoint(time=time, by_population=by_population)

    def simulate(self, tidy: bool = False):
        """Run the simulation and return sampled multi-population state points or a tidy DataFrame if tidy=True."""

        dt = self.simulation.time_step
        steps = int(round(self.simulation.duration / dt))
        state = self._initial_state()

        trajectory = [self._structured_state_point(0.0, state)]

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

            raw_state = self.stepper.step(state=state, time_step=dt, derivatives=self._derivatives)
            state = self._normalize_state(raw_state)
            if step % self.simulation.output_stride == 0:
                trajectory.append(self._structured_state_point(step * dt, state))

        if not tidy:
            return trajectory

        # Tidy DataFrame: one row per time, population, compartment
        import pandas as pd
        records = []
        for pt in trajectory:
            for pop_name, pop_state in pt.by_population.items():
                records.append({
                    "time": pt.time,
                    "population": pop_name,
                    "susceptible": pop_state.susceptible,
                    "exposed": pop_state.exposed,
                    "infected": pop_state.infected,
                    "recovered": pop_state.recovered,
                    "total_population": pop_state.total_population,
                })
        return pd.DataFrame(records)

    def contact_matrix_df(self) -> pd.DataFrame:
        """Return the contact matrix as a pandas DataFrame with population names as index/columns."""
        return pd.DataFrame(self.contact_matrix, index=self._population_names, columns=self._population_names)
