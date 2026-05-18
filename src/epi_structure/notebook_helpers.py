"""Helpers to keep pedagogy notebooks focused on narrative over mechanics."""

from __future__ import annotations

import importlib
import sys
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle, FancyArrowPatch

from .model import StructuredEpidemicModel, StructuredStatePoint, TransitionProbabilityStepper
from .scenarios import Scenario, get_scenario


def reload_package_modules(package_name: str = "epi_structure"):
    """Reload a package and all currently loaded submodules under it.

    This is primarily intended for notebook workflows where iterative edits
    should be reflected without maintaining a manual list of module reloads.
    """

    root = importlib.import_module(package_name)
    module_names = [
        name
        for name in sys.modules
        if name == package_name or name.startswith(f"{package_name}.")
    ]
    module_names.sort(key=lambda name: name.count("."), reverse=True)

    for module_name in module_names:
        module = sys.modules.get(module_name)
        if module is None:
            continue
        importlib.reload(module)

    return importlib.import_module(package_name)


def run_scenario(
    scenario: str | Scenario,
    simulation=None,
    contact_matrix_override: list[list[float]] | None = None,
    intervention_plan=None,
    engine: str = "rk4",
    seed: int | None = None,
) -> list[StructuredStatePoint]:
    """Run a scenario by name or object, with optional simulation and matrix overrides.

    Use ``engine="probability"`` to sample whole-number transitions instead of
    integrating fractional compartment counts.
    """
    scenario_obj = get_scenario(scenario) if isinstance(scenario, str) else scenario
    if not isinstance(scenario_obj, Scenario) and hasattr(scenario_obj, "build"):
        scenario_obj = scenario_obj.build()

    stepper = None
    normalized_engine = engine.strip().lower()
    if normalized_engine in {"rk4", "deterministic"}:
        stepper = None
    elif normalized_engine in {"probability", "transition_probability", "stochastic"}:
        stepper = TransitionProbabilityStepper(seed=seed)
    else:
        raise ValueError("engine must be 'rk4' or 'probability'")

    model = StructuredEpidemicModel(
        populations=scenario_obj.populations,
        contact_matrix=contact_matrix_override or scenario_obj.contact_matrix,
        simulation=simulation or scenario_obj.simulation,
        intervention_plan=(scenario_obj.intervention_plan if intervention_plan is None else intervention_plan),
        stepper=stepper,
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


def trajectory_to_df(trajectory: list[StructuredStatePoint]) -> pd.DataFrame:
    """Convert structured trajectory output to a tidy DataFrame with incidence proxy."""
    records: list[dict[str, float | str]] = []
    for point in trajectory:
        for pop_name, state in point.by_population.items():
            records.append(
                {
                    "time": point.time,
                    "population": pop_name,
                    "S": state.susceptible,
                    "E": state.exposed,
                    "I": state.infected,
                    "R": state.recovered,
                    "N": state.total_population,
                }
            )

    df = pd.DataFrame(records).sort_values(["population", "time"]).reset_index(drop=True)
    df["new_exposed_interval"] = (-df.groupby("population")["S"].diff().fillna(0.0)).clip(lower=0.0)
    return df


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


def show_params(name: str, **kwargs) -> None:
    """Print a compact parameter summary block for notebooks."""
    print(f"\n{name}:")
    for key, value in kwargs.items():
        if isinstance(value, float):
            print(f"  {key:24s} = {value:.6f}")
        else:
            print(f"  {key:24s} = {value}")


def simulate_sir_euler(
    N: int,
    I0: float,
    R0_value: float,
    gamma: float,
    duration: float,
    dt: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Forward-Euler SIR simulation returning times, S, I, R, incidence."""
    beta = R0_value * gamma
    S0 = N - I0
    R = 0.0

    S_traj = [float(S0)]
    I_traj = [float(I0)]
    R_traj = [R]
    incidence: list[float] = []

    S = float(S0)
    I = float(I0)

    for _ in range(int(duration / dt)):
        force_of_infection = beta * S * I / N
        dS = -force_of_infection
        dI = force_of_infection - gamma * I
        dR = gamma * I

        S += dS * dt
        I += dI * dt
        R += dR * dt

        S_traj.append(max(0.0, S))
        I_traj.append(max(0.0, I))
        R_traj.append(max(0.0, R))
        incidence.append(max(0.0, force_of_infection * dt))

    times = np.arange(0, duration + dt, dt)
    return (
        times[: len(S_traj)],
        np.array(S_traj),
        np.array(I_traj),
        np.array(R_traj),
        np.array(incidence),
    )


def simulate_seir_euler(
    N: int,
    I0: float,
    E0: float,
    R0_value: float,
    gamma: float,
    sigma: float,
    duration: float,
    dt: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Forward-Euler SEIR simulation returning times, S, E, I, R, incidence."""
    beta = R0_value * gamma
    S0 = N - I0 - E0
    R = 0.0

    S_traj = [float(S0)]
    E_traj = [float(E0)]
    I_traj = [float(I0)]
    R_traj = [R]
    incidence: list[float] = []

    S = float(S0)
    E = float(E0)
    I = float(I0)

    for _ in range(int(duration / dt)):
        force_of_infection = beta * S * I / N
        dS = -force_of_infection
        dE = force_of_infection - sigma * E
        dI = sigma * E - gamma * I
        dR = gamma * I

        S += dS * dt
        E += dE * dt
        I += dI * dt
        R += dR * dt

        S_traj.append(max(0.0, S))
        E_traj.append(max(0.0, E))
        I_traj.append(max(0.0, I))
        R_traj.append(max(0.0, R))
        incidence.append(max(0.0, force_of_infection * dt))

    times = np.arange(0, duration + dt, dt)
    return (
        times[: len(S_traj)],
        np.array(S_traj),
        np.array(E_traj),
        np.array(I_traj),
        np.array(R_traj),
        np.array(incidence),
    )


def simulate_seirs_euler(
    N: int,
    I0: float,
    E0: float,
    R0_value: float,
    gamma: float,
    sigma: float,
    rho: float,
    duration: float,
    dt: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Forward-Euler SEIRS simulation returning times, S, E, I, R, incidence."""
    beta = R0_value * gamma

    nsteps = int(duration / dt)
    times = np.arange(nsteps) * dt

    S = np.zeros(nsteps)
    E = np.zeros(nsteps)
    I = np.zeros(nsteps)
    R = np.zeros(nsteps)
    incidence = np.zeros(nsteps)

    S[0] = N - I0 - E0
    E[0] = E0
    I[0] = I0
    R[0] = 0.0

    for step in range(1, nsteps):
        S_curr, E_curr, I_curr, R_curr = S[step - 1], E[step - 1], I[step - 1], R[step - 1]

        dS_dt = -beta * S_curr * I_curr / N + rho * R_curr
        dE_dt = beta * S_curr * I_curr / N - sigma * E_curr
        dI_dt = sigma * E_curr - gamma * I_curr
        dR_dt = gamma * I_curr - rho * R_curr

        S[step] = S_curr + dS_dt * dt
        E[step] = E_curr + dE_dt * dt
        I[step] = I_curr + dI_dt * dt
        R[step] = R_curr + dR_dt * dt

        incidence[step] = sigma * E_curr * dt

    return times, S, E, I, R, incidence


def plot_sir_trajectory(
    times: np.ndarray,
    S: np.ndarray,
    I: np.ndarray,
    R: np.ndarray,
    incidence: np.ndarray,
    title: str = "",
):
    """Plot SIR trajectories and incidence in a 1x2 panel."""
    fig = plt.figure(figsize=(14, 5))
    gs = GridSpec(1, 2, figure=fig)

    ax0 = fig.add_subplot(gs[0, 0])
    ax0.plot(times, S, label="S (Susceptible)", linewidth=2, color="#1f77b4")
    ax0.plot(times, I, label="I (Infected)", linewidth=2, color="#ff7f0e")
    ax0.plot(times, R, label="R (Recovered)", linewidth=2, color="#2ca02c")
    ax0.set_xlabel("Time (days)", fontsize=11)
    ax0.set_ylabel("Count", fontsize=11)
    ax0.set_title("Compartment Trajectories", fontsize=12, fontweight="bold")
    ax0.legend(loc="best")
    ax0.grid(True, alpha=0.3)

    ax1 = fig.add_subplot(gs[0, 1])
    ax1.bar(times[:-1], incidence, width=1, alpha=0.7, color="#d62728", label="Daily incidence")
    ax1.set_xlabel("Time (days)", fontsize=11)
    ax1.set_ylabel("New infections", fontsize=11)
    ax1.set_title("Incidence (New Infections per Day)", fontsize=12, fontweight="bold")
    ax1.legend(loc="best")
    ax1.grid(True, alpha=0.3, axis="y")

    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold")
    plt.tight_layout()
    return fig


def plot_seir_trajectory(
    times: np.ndarray,
    S: np.ndarray,
    E: np.ndarray,
    I: np.ndarray,
    R: np.ndarray,
    incidence: np.ndarray,
    title: str = "",
):
    """Plot SEIR trajectories and incidence in a 1x2 panel."""
    fig = plt.figure(figsize=(14, 5))
    gs = GridSpec(1, 2, figure=fig)

    ax0 = fig.add_subplot(gs[0, 0])
    ax0.plot(times, S, label="S (Susceptible)", linewidth=2, color="#1f77b4")
    ax0.plot(times, E, label="E (Exposed)", linewidth=2, color="#ff7f0e")
    ax0.plot(times, I, label="I (Infected)", linewidth=2, color="#d62728")
    ax0.plot(times, R, label="R (Recovered)", linewidth=2, color="#2ca02c")
    ax0.set_xlabel("Time (days)", fontsize=11)
    ax0.set_ylabel("Count", fontsize=11)
    ax0.set_title("Compartment Trajectories", fontsize=12, fontweight="bold")
    ax0.legend(loc="best", fontsize=10)
    ax0.grid(True, alpha=0.3)

    ax1 = fig.add_subplot(gs[0, 1])
    ax1.bar(times[:-1], incidence, width=1, alpha=0.7, color="#ff7f0e", label="Daily incidence")
    ax1.set_xlabel("Time (days)", fontsize=11)
    ax1.set_ylabel("New infections", fontsize=11)
    ax1.set_title("Incidence (New Infections per Day)", fontsize=12, fontweight="bold")
    ax1.legend(loc="best")
    ax1.grid(True, alpha=0.3, axis="y")

    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold")
    plt.tight_layout()
    return fig


def plot_seirs_trajectory(
    times: np.ndarray,
    S: np.ndarray,
    E: np.ndarray,
    I: np.ndarray,
    R: np.ndarray,
    incidence: np.ndarray,
    title: str = "",
):
    """Plot SEIRS trajectories and incidence in a 1x2 panel."""
    fig = plt.figure(figsize=(14, 5))
    gs = GridSpec(1, 2, figure=fig)

    ax0 = fig.add_subplot(gs[0, 0])
    ax0.plot(times, S, label="S (Susceptible)", linewidth=2, color="#1f77b4")
    ax0.plot(times, E, label="E (Exposed)", linewidth=2, color="#ff7f0e")
    ax0.plot(times, I, label="I (Infected)", linewidth=2, color="#d62728")
    ax0.plot(times, R, label="R (Recovered)", linewidth=2, color="#2ca02c")
    ax0.set_xlabel("Time (days)", fontsize=11)
    ax0.set_ylabel("Number of Individuals", fontsize=11)
    ax0.legend(loc="best", fontsize=10)
    ax0.set_title("Compartment Trajectories", fontsize=12, fontweight="bold")
    ax0.grid(True, alpha=0.3)

    ax1 = fig.add_subplot(gs[0, 1])
    ax1.bar(times, incidence, width=0.8, color="#d62728", alpha=0.7, edgecolor="black", linewidth=0.5)
    ax1.set_xlabel("Time (days)", fontsize=11)
    ax1.set_ylabel("Daily Incidence", fontsize=11)
    ax1.set_title("Daily New Infections (E→I flow)", fontsize=12, fontweight="bold")
    ax1.grid(True, alpha=0.3, axis="y")

    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.show()
    return fig
