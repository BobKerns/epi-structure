"""Helpers to keep pedagogy notebooks focused on narrative over mechanics."""

from __future__ import annotations

from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Circle, FancyArrowPatch

from .model import StructuredEpidemicModel, StructuredStatePoint
from .scenarios import get_scenario


def run_scenario(scenario_name: str) -> list[StructuredStatePoint]:
    """Load and run a named scenario."""
    scenario = get_scenario(scenario_name)
    model = StructuredEpidemicModel(
        populations=scenario.populations,
        contact_matrix=scenario.contact_matrix,
        simulation=scenario.simulation,
        intervention_plan=scenario.intervention_plan,
    )
    return model.simulate()


def get_intervention_plan(scenario_name: str):
    """Get the intervention plan for a scenario."""
    return get_scenario(scenario_name).intervention_plan


def get_contact_matrix(scenario_name: str) -> list[list[float]]:
    """Get the base contact matrix for a scenario."""
    return get_scenario(scenario_name).contact_matrix


def effective_matrix(
    scenario_name: str,
    base_contact_matrix: list[list[float]],
    population_names: list[str],
    time: float,
) -> list[list[float]]:
    """Compute effective contact matrix at a given time under a scenario."""
    plan = get_intervention_plan(scenario_name)
    if plan is None:
        return [row[:] for row in base_contact_matrix]
    return plan.matrix_at_time(base_contact_matrix, population_names, time)


def intervention_targets(
    scenario_name: str,
    population_names: list[str],
) -> tuple[set[str], set[tuple[str, str]]]:
    """Extract targeted groups and bridge links from a scenario's intervention plan."""
    plan = get_intervention_plan(scenario_name)
    groups: set[str] = set()
    links: set[tuple[str, str]] = set()

    if plan is None:
        return groups, links

    for intervention in plan.interventions:
        if isinstance(intervention.within_scale, dict):
            groups.update(intervention.within_scale.keys())
        elif isinstance(intervention.within_scale, (int, float)):
            groups.update(population_names)

        if isinstance(intervention.bridge_scale, dict):
            links.update(intervention.bridge_scale.keys())
        elif isinstance(intervention.bridge_scale, (int, float)):
            for to_name in population_names:
                for from_name in population_names:
                    if to_name != from_name:
                        links.add((to_name, from_name))

    return groups, links


def _fmt_weight(value: float) -> str:
    abs_value = abs(value)
    if abs_value == 0:
        return "0"
    if abs_value < 0.01:
        return f"{value:.2e}"
    return f"{value:.3f}"


def draw_group_bridge_map(
    ax,
    names: list[str],
    matrix: list[list[float]],
    title: str,
    highlighted_groups: set[str] | None = None,
    highlighted_links: set[tuple[str, str]] | None = None,
) -> None:
    """Draw group nodes and directed bridge links from a contact matrix."""
    highlighted_groups = highlighted_groups or set()
    highlighted_links = highlighted_links or set()
    index = {name: i for i, name in enumerate(names)}
    positions = {name: (float(i), 0.0) for i, name in enumerate(names)}

    max_entry = max(max(row) for row in matrix) if matrix else 1.0
    max_entry = max(max_entry, 1e-9)

    for name, (x, y) in positions.items():
        i = index[name]
        diag = matrix[i][i]
        radius = 0.16 + 0.08 * (diag / max_entry)
        edgecolor = "tomato" if name in highlighted_groups else "black"
        circle = Circle((x, y), radius=radius, facecolor="#d9edf7", edgecolor=edgecolor, linewidth=2)
        ax.add_patch(circle)
        ax.text(
            x,
            y + radius + 0.08,
            f"{name}\nwithin={_fmt_weight(diag)}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    n = len(names)
    for to_idx in range(n):
        for from_idx in range(n):
            if to_idx == from_idx:
                continue

            to_name = names[to_idx]
            from_name = names[from_idx]
            value = matrix[to_idx][from_idx]
            start = positions[from_name]
            end = positions[to_name]

            if start == end:
                continue

            is_highlighted = (to_name, from_name) in highlighted_links
            color = "tomato" if is_highlighted else "#777777"
            linewidth = 1.0 + 3.0 * (value / max_entry)

            arrow = FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=12,
                color=color,
                linewidth=linewidth,
                connectionstyle="arc3,rad=0.2" if from_idx < to_idx else "arc3,rad=-0.2",
                shrinkA=20,
                shrinkB=20,
            )
            ax.add_patch(arrow)

            mid_x = (start[0] + end[0]) / 2
            mid_y = (start[1] + end[1]) / 2
            offset = 0.16 if from_idx < to_idx else -0.16
            ax.text(
                mid_x,
                mid_y + offset,
                _fmt_weight(value),
                ha="center",
                va="center",
                fontsize=8,
                color=color,
            )

    ax.set_title(title)
    ax.set_xlim(-0.6, len(names) - 0.4)
    ax.set_ylim(-0.8, 0.8)
    ax.set_aspect("equal")
    ax.axis("off")


def plot_matrix_transformation(
    names: list[str],
    base_matrix: list[list[float]],
    effective: list[list[float]],
    title: str,
    highlighted_groups: set[str] | None = None,
    highlighted_links: set[tuple[str, str]] | None = None,
):
    """Plot base and effective contact maps side-by-side."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    draw_group_bridge_map(
        axes[0],
        names,
        base_matrix,
        "Base Contact Map",
    )
    draw_group_bridge_map(
        axes[1],
        names,
        effective,
        title,
        highlighted_groups=highlighted_groups,
        highlighted_links=highlighted_links,
    )
    fig.suptitle("Selected Subset -> Matrix Transformation")
    plt.tight_layout()
    return fig, axes


def summarize_two_population_runs(
    scenario_names: Iterable[str],
    runs: dict[str, list[StructuredStatePoint]],
    general_name: str = "general",
    cluster_name: str = "cluster",
    general_size: float = 10_000,
    cluster_size: float = 500,
) -> pd.DataFrame:
    """Summarize peak and endpoint metrics for two named populations."""

    def summarize_run(name: str, trajectory: list[StructuredStatePoint]) -> dict[str, float | str]:
        general = [pt.by_population[general_name].infected for pt in trajectory]
        cluster = [pt.by_population[cluster_name].infected for pt in trajectory]
        times = [pt.time for pt in trajectory]

        general_pct = [100.0 * x / general_size for x in general]
        cluster_pct = [100.0 * x / cluster_size for x in cluster]

        g_peak_idx = max(range(len(general)), key=lambda i: general[i])
        c_peak_idx = max(range(len(cluster)), key=lambda i: cluster[i])

        return {
            "scenario": name,
            "general_peak_I": general[g_peak_idx],
            "general_peak_t": times[g_peak_idx],
            "general_end_I": general[-1],
            "general_peak_pct": general_pct[g_peak_idx],
            "cluster_peak_I": cluster[c_peak_idx],
            "cluster_peak_t": times[c_peak_idx],
            "cluster_end_I": cluster[-1],
            "cluster_peak_pct": cluster_pct[c_peak_idx],
        }

    return pd.DataFrame([summarize_run(name, runs[name]) for name in scenario_names])


def plot_population_infected_comparison(
    runs: dict[str, list[StructuredStatePoint]],
    population_name: str,
    population_size: float,
    suptitle: str,
):
    """Plot infected count and percent trajectories for one population across scenarios."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=True)

    for name, trajectory in runs.items():
        times = [pt.time for pt in trajectory]
        infected = [pt.by_population[population_name].infected for pt in trajectory]
        infected_pct = [100.0 * x / population_size for x in infected]

        axes[0].plot(times, infected, label=name)
        axes[1].plot(times, infected_pct, label=name)

    axes[0].set_title(f"{population_name.title()} Population: Infected Count")
    axes[0].set_xlabel("time")
    axes[0].set_ylabel("infected count")
    axes[0].grid(True)
    axes[0].legend()

    axes[1].set_title(f"{population_name.title()} Population: Percent Infected")
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("infected (%)")
    axes[1].grid(True)
    axes[1].legend()

    fig.suptitle(suptitle)
    plt.tight_layout()
    return fig, axes
