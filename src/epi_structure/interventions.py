"""Intervention helpers for structured epidemic simulations."""

from dataclasses import dataclass


BridgeLink = tuple[str, str]


@dataclass(frozen=True, slots=True)
class ContactMatrixIntervention:
    """Time-bounded scaling intervention applied to contact matrix entries.

    Attributes:
        start_time: Simulation time (inclusive) when intervention becomes active.
        end_time: Simulation time (exclusive) when intervention stops; None means no end.
        within_scale: Scalar or per-population scale for diagonal entries.
        bridge_scale: Scalar or per-(to, from) scale for off-diagonal entries.
    """

    start_time: float
    end_time: float | None = None
    within_scale: float | dict[str, float] | None = None
    bridge_scale: float | dict[tuple[str, str], float] | None = None

    def __post_init__(self) -> None:
        if self.start_time < 0:
            raise ValueError("start_time must be non-negative")
        if self.end_time is not None and self.end_time <= self.start_time:
            raise ValueError("end_time must be greater than start_time")

        for value in _iter_scales(self.within_scale):
            if value < 0:
                raise ValueError("within_scale values must be non-negative")
        for value in _iter_scales(self.bridge_scale):
            if value < 0:
                raise ValueError("bridge_scale values must be non-negative")

    def is_active(self, time: float) -> bool:
        """Return whether this intervention is active at a simulation time."""

        if time < self.start_time:
            return False
        if self.end_time is None:
            return True
        return time < self.end_time


def _iter_scales(scale: float | dict[object, float] | None) -> list[float]:
    if scale is None:
        return []
    if isinstance(scale, (int, float)):
        return [float(scale)]
    return [float(value) for value in scale.values()]


class InterventionPlan:
    """Collection of interventions that transforms a base contact matrix over time."""

    def __init__(self, interventions: list[ContactMatrixIntervention] | None = None) -> None:
        self.interventions = sorted(interventions or [], key=lambda item: item.start_time)

    def matrix_at_time(
        self,
        base_contact_matrix: list[list[float]],
        population_names: list[str],
        time: float,
    ) -> list[list[float]]:
        """Return effective contact matrix at time after applying active interventions."""

        matrix = [row[:] for row in base_contact_matrix]
        index = {name: idx for idx, name in enumerate(population_names)}

        for intervention in self.interventions:
            if not intervention.is_active(time):
                continue

            self._apply_within_scale(matrix, index, intervention.within_scale)
            self._apply_bridge_scale(matrix, index, intervention.bridge_scale)

        return matrix

    @staticmethod
    def _apply_within_scale(
        matrix: list[list[float]],
        index: dict[str, int],
        within_scale: float | dict[str, float] | None,
    ) -> None:
        if within_scale is None:
            return

        if isinstance(within_scale, (int, float)):
            scale = float(within_scale)
            for i in range(len(matrix)):
                matrix[i][i] *= scale
            return

        for name, scale in within_scale.items():
            if name not in index:
                raise ValueError(f"Unknown population in within_scale: {name}")
            i = index[name]
            matrix[i][i] *= float(scale)

    @staticmethod
    def _apply_bridge_scale(
        matrix: list[list[float]],
        index: dict[str, int],
        bridge_scale: float | dict[tuple[str, str], float] | None,
    ) -> None:
        if bridge_scale is None:
            return

        if isinstance(bridge_scale, (int, float)):
            scale = float(bridge_scale)
            for i in range(len(matrix)):
                for j in range(len(matrix)):
                    if i != j:
                        matrix[i][j] *= scale
            return

        for (to_group, from_group), scale in bridge_scale.items():
            if to_group not in index or from_group not in index:
                raise ValueError(f"Unknown population in bridge_scale: {(to_group, from_group)}")
            i = index[to_group]
            j = index[from_group]
            if i == j:
                raise ValueError("bridge_scale entries must be off-diagonal")
            matrix[i][j] *= float(scale)


def intervene_within_groups(
    groups: list[str],
    scale: float,
    *,
    start_time: float,
    end_time: float | None = None,
) -> ContactMatrixIntervention:
    """Create a within-group intervention for a selected set of groups.

    This follows a simple functional pattern: select model subset (groups), then
    apply a transformation (diagonal scale) over a time window.
    """

    if not groups:
        raise ValueError("groups must not be empty")

    return ContactMatrixIntervention(
        start_time=start_time,
        end_time=end_time,
        within_scale={group: float(scale) for group in groups},
    )


def intervene_bridge_links(
    links: list[BridgeLink],
    scale: float,
    *,
    start_time: float,
    end_time: float | None = None,
    symmetric: bool = False,
) -> ContactMatrixIntervention:
    """Create a bridge-link intervention for selected directed links.

    If ``symmetric`` is True, each selected link also applies to its reverse.
    """

    if not links:
        raise ValueError("links must not be empty")

    bridge_scale: dict[BridgeLink, float] = {}
    for to_group, from_group in links:
        if to_group == from_group:
            raise ValueError("bridge links must be off-diagonal")
        bridge_scale[(to_group, from_group)] = float(scale)
        if symmetric:
            bridge_scale[(from_group, to_group)] = float(scale)

    return ContactMatrixIntervention(
        start_time=start_time,
        end_time=end_time,
        bridge_scale=bridge_scale,
    )


def compose_intervention_plan(*interventions: ContactMatrixIntervention) -> InterventionPlan:
    """Compose an intervention plan from one or more intervention transforms."""

    return InterventionPlan(list(interventions))
