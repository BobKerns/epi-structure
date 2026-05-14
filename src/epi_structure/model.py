"""Core epidemic model implementation."""

from dataclasses import dataclass
from typing import Callable, Protocol

from .parameters import PopulationParameters, SimulationParameters


StateVector = tuple[float, float, float, float]
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
        return (
            base[0] + delta[0] * scale,
            base[1] + delta[1] * scale,
            base[2] + delta[2] * scale,
            base[3] + delta[3] * scale,
        )

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

        return (
            state[0] + time_step * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0,
            state[1] + time_step * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0,
            state[2] + time_step * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]) / 6.0,
            state[3] + time_step * (k1[3] + 2.0 * k2[3] + 2.0 * k3[3] + k4[3]) / 6.0,
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

    def simulate(self) -> list[StatePoint]:
        """Run the simulation and return sampled state points."""

        dt = self.simulation.time_step
        steps = int(round(self.simulation.duration / dt))
        state: StateVector = tuple(float(v) for v in self.population.initial_state)

        trajectory = [StatePoint(0.0, *state)]

        for step in range(1, steps + 1):
            raw_state = self.stepper.step(state=state, time_step=dt, derivatives=self._derivatives)
            state = self._normalize_state(raw_state, float(self.population.size))
            if step % self.simulation.output_stride == 0:
                trajectory.append(StatePoint(step * dt, *state))

        return trajectory
