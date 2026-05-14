"""Helpers for constructing contact matrices from lightweight named inputs."""


def build_contact_matrix(
    pop_names: list[str],
    within: float | dict[str, float],
    bridges: dict[tuple[str, str], float] | None = None,
) -> list[list[float]]:
    """Build a square contact matrix from named inputs.

    Matrix meaning: beta_ij is transmission contribution from group j to group i.

    Args:
        pop_names: Population names in the simulation order expected by the model.
        within: Either one scalar for all diagonal entries or a dict of per-population diagonals.
        bridges: Sparse off-diagonal overrides keyed by (to_group, from_group).

    Returns:
        A square matrix with shape (n_populations, n_populations).
    """

    if not pop_names:
        raise ValueError("pop_names must not be empty")

    bridges = bridges or {}
    index = {name: i for i, name in enumerate(pop_names)}
    n = len(pop_names)
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]

    if isinstance(within, (int, float)):
        if within < 0:
            raise ValueError("within must be non-negative")
        for i in range(n):
            matrix[i][i] = float(within)
    else:
        for name, value in within.items():
            if name not in index:
                raise ValueError(f"Unknown population in within: {name}")
            if value < 0:
                raise ValueError("within values must be non-negative")
            matrix[index[name]][index[name]] = float(value)

    for (to_group, from_group), value in bridges.items():
        if to_group not in index or from_group not in index:
            raise ValueError(f"Unknown population in bridge: {(to_group, from_group)}")
        if value < 0:
            raise ValueError("bridge values must be non-negative")
        matrix[index[to_group]][index[from_group]] = float(value)

    return matrix


def preset_general_cluster() -> tuple[list[str], list[list[float]]]:
    """Return a default two-population general/cluster contact matrix preset."""

    names = ["general", "cluster"]
    within = {"general": 0.00004, "cluster": 0.0005}
    bridges = {
        ("general", "cluster"): 0.00018,
        ("cluster", "general"): 0.00001,
    }
    return names, build_contact_matrix(names, within=within, bridges=bridges)
