# Epidemic Population Structure Demo

Interactive demonstration of how population structure drives epidemic dynamics in ways homogeneous models miss.

Current focus: Phase 1 notebook-first MVP for model development, explanation, and validation.

## Why This Exists

The project demonstrates a specific claim with policy relevance:

- A general population with R0 < 1 does not sustain transmission on its own.
- A supercritical cluster can sustain local transmission.
- A bridge from that cluster can continuously reseed the general population.
- Targeted interventions on clusters and bridges can clear general-population
  incidence without changing general-population R0.

For design rationale and narrative context, see DESIGN.md.

## Project Status

Early build phase.

Implemented now:
- Design documentation
- Source artifacts (thread and whiteboard image)
- License

Planned next:
- Python package scaffold
- Notebook sequence for stages 1 to 6
- Baseline parameter table and validation scenarios

## Planned Repository Layout

```text
src/epi_structure/
  __init__.py
  model.py
  parameters.py
  interventions.py
notebooks/
  01_sir_single.ipynb
  02_seir.ipynb
  03_seirs.ipynb
  04_two_populations.ipynb
  05_interventions.ipynb
DESIGN.md
README.md
LICENSE
```

## Model Scope (Phase 1)

- SEIRS dynamics with a contact matrix for between-group transmission
- Density-dependent transmission form consistent with Kermack-McKendrick equations
- Deterministic RK4 time integration
- Intervention toggles for within-cluster and bridge transmission scaling

Out of scope for Phase 1:
- Public Observable deployment
- Pyodide browser execution
- Full multi-cluster stage 7 model

## Quick Start (Will Be Updated As Code Lands)

### 1) Create environment

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
```

### 2) Install dependencies

```bash
pip install -e .
```

### 3) Run notebooks

```bash
jupyter lab
```

Note: these commands are placeholders until pyproject.toml and notebooks are added.

## Data and Narrative Artifacts

- thread.md: source thread motivating the project
- images/whiteboard.png: original whiteboard SIR equations photo

## References

Primary references are maintained in DESIGN.md to keep design and citation context together during Phase 1.

## License

MIT. See LICENSE.
