# Notebooks

This directory contains Jupyter notebooks for exploring and demonstrating the epidemic population structure model.

## Current Notebooks

### Pedagogical Sequence (Foundation-First)

**01_sir.ipynb**: SIR model with R₀ sensitivity analysis
- Simplest compartmental model (no latent period, permanent immunity)
- Baseline dynamics and visualization
- Effective reproduction number R_t over time
- **Teaching focus**: R₀ sensitivity, intuition failures (nonlinear outcomes), extrapolation failures
- Demonstrates how small changes in R₀ produce large changes in peak infections and attack rate
- Early exponential fit vs full trajectory comparison

**02_seir.ipynb**: SEIR model adds latent period
- S → E → I → R progression
- How latent period affects timing and peak while R₀ remains constant
- Comparison with SIR baseline

**03_seirs.ipynb**: SEIRS adds waning immunity
- Full loop with loss of immunity
- Endemic equilibrium and cycle behavior
- How R₀ relates to endemic prevalence

**04_seird.ipynb**: SEIRD model adds mortality
- S → E → I → R plus disease-induced deaths
- CFR-governed mortality flow and conservation checks
- Comparison with SEIR timing and severity outcomes

**05_two_population.ipynb**: Multi-population dynamics
- Two populations (general + cluster) with contact matrix
- Within-group and bridge transmission
- How coupling affects each population differently
- Effective R₀ is population-dependent

**06_superspreader_party_crawl.ipynb**: Sparse coupling and re-seeding dynamics
- General-only subcritical regime ($R_0 < 1$) contrasted with coupled cluster dynamics
- Superspreader cluster surges and spillback into the general population

**07_rt_wells_riley_decomposition.ipynb**: Mechanistic $R_t$ decomposition
- Decomposes $R_t$ into contact rate, infectious duration, susceptibility, and per-contact infection probability
- Uses a Wells-Riley style dose model for infection probability terms

**08_interventions.ipynb**: Policy language to matrix-scaling transformations
- Ventilation, tracing, distancing as transmission scalar interventions
- Group/bridge link visualization with subset-selection + transformation pattern
- Scenario comparison (baseline, cluster ventilation, bridge tracing, targeted combo)

**09_timed_interventions.ipynb**: Time-varying interventions
- Component-level policy schedules for $c_t$ and Wells-Riley terms
- Time-bounded matrix interventions in the structured model

## Public Access

If GitHub fails to render a notebook, use nbviewer. For cloud execution without local
setup, use Colab or Binder.

| Notebook | GitHub | nbviewer | Colab | Binder |
| --- | --- | --- | --- | --- |
| 01_sir.ipynb | [view](https://github.com/BobKerns/epi-structure/blob/main/notebooks/01_sir.ipynb) | [view](https://nbviewer.org/github/BobKerns/epi-structure/blob/main/notebooks/01_sir.ipynb) | [run](https://colab.research.google.com/github/BobKerns/epi-structure/blob/main/notebooks/01_sir.ipynb) | [run](https://mybinder.org/v2/gh/BobKerns/epi-structure/main?labpath=notebooks/01_sir.ipynb) |
| 02_seir.ipynb | [view](https://github.com/BobKerns/epi-structure/blob/main/notebooks/02_seir.ipynb) | [view](https://nbviewer.org/github/BobKerns/epi-structure/blob/main/notebooks/02_seir.ipynb) | [run](https://colab.research.google.com/github/BobKerns/epi-structure/blob/main/notebooks/02_seir.ipynb) | [run](https://mybinder.org/v2/gh/BobKerns/epi-structure/main?labpath=notebooks/02_seir.ipynb) |
| 03_seirs.ipynb | [view](https://github.com/BobKerns/epi-structure/blob/main/notebooks/03_seirs.ipynb) | [view](https://nbviewer.org/github/BobKerns/epi-structure/blob/main/notebooks/03_seirs.ipynb) | [run](https://colab.research.google.com/github/BobKerns/epi-structure/blob/main/notebooks/03_seirs.ipynb) | [run](https://mybinder.org/v2/gh/BobKerns/epi-structure/main?labpath=notebooks/03_seirs.ipynb) |
| 04_seird.ipynb | [view](https://github.com/BobKerns/epi-structure/blob/main/notebooks/04_seird.ipynb) | [view](https://nbviewer.org/github/BobKerns/epi-structure/blob/main/notebooks/04_seird.ipynb) | [run](https://colab.research.google.com/github/BobKerns/epi-structure/blob/main/notebooks/04_seird.ipynb) | [run](https://mybinder.org/v2/gh/BobKerns/epi-structure/main?labpath=notebooks/04_seird.ipynb) |
| 05_two_population.ipynb | [view](https://github.com/BobKerns/epi-structure/blob/main/notebooks/05_two_population.ipynb) | [view](https://nbviewer.org/github/BobKerns/epi-structure/blob/main/notebooks/05_two_population.ipynb) | [run](https://colab.research.google.com/github/BobKerns/epi-structure/blob/main/notebooks/05_two_population.ipynb) | [run](https://mybinder.org/v2/gh/BobKerns/epi-structure/main?labpath=notebooks/05_two_population.ipynb) |
| 06_superspreader_party_crawl.ipynb | [view](https://github.com/BobKerns/epi-structure/blob/main/notebooks/06_superspreader_party_crawl.ipynb) | [view](https://nbviewer.org/github/BobKerns/epi-structure/blob/main/notebooks/06_superspreader_party_crawl.ipynb) | [run](https://colab.research.google.com/github/BobKerns/epi-structure/blob/main/notebooks/06_superspreader_party_crawl.ipynb) | [run](https://mybinder.org/v2/gh/BobKerns/epi-structure/main?labpath=notebooks/06_superspreader_party_crawl.ipynb) |
| 07_rt_wells_riley_decomposition.ipynb | [view](https://github.com/BobKerns/epi-structure/blob/main/notebooks/07_rt_wells_riley_decomposition.ipynb) | [view](https://nbviewer.org/github/BobKerns/epi-structure/blob/main/notebooks/07_rt_wells_riley_decomposition.ipynb) | [run](https://colab.research.google.com/github/BobKerns/epi-structure/blob/main/notebooks/07_rt_wells_riley_decomposition.ipynb) | [run](https://mybinder.org/v2/gh/BobKerns/epi-structure/main?labpath=notebooks/07_rt_wells_riley_decomposition.ipynb) |
| 08_interventions.ipynb | [view](https://github.com/BobKerns/epi-structure/blob/main/notebooks/08_interventions.ipynb) | [view](https://nbviewer.org/github/BobKerns/epi-structure/blob/main/notebooks/08_interventions.ipynb) | [run](https://colab.research.google.com/github/BobKerns/epi-structure/blob/main/notebooks/08_interventions.ipynb) | [run](https://mybinder.org/v2/gh/BobKerns/epi-structure/main?labpath=notebooks/08_interventions.ipynb) |
| 09_timed_interventions.ipynb | [view](https://github.com/BobKerns/epi-structure/blob/main/notebooks/09_timed_interventions.ipynb) | [view](https://nbviewer.org/github/BobKerns/epi-structure/blob/main/notebooks/09_timed_interventions.ipynb) | [run](https://colab.research.google.com/github/BobKerns/epi-structure/blob/main/notebooks/09_timed_interventions.ipynb) | [run](https://mybinder.org/v2/gh/BobKerns/epi-structure/main?labpath=notebooks/09_timed_interventions.ipynb) |
| contact_matrix_exploration.ipynb | [view](https://github.com/BobKerns/epi-structure/blob/main/notebooks/contact_matrix_exploration.ipynb) | [view](https://nbviewer.org/github/BobKerns/epi-structure/blob/main/notebooks/contact_matrix_exploration.ipynb) | [run](https://colab.research.google.com/github/BobKerns/epi-structure/blob/main/notebooks/contact_matrix_exploration.ipynb) | [run](https://mybinder.org/v2/gh/BobKerns/epi-structure/main?labpath=notebooks/contact_matrix_exploration.ipynb) |
| exploration.ipynb | [view](https://github.com/BobKerns/epi-structure/blob/main/notebooks/exploration.ipynb) | [view](https://nbviewer.org/github/BobKerns/epi-structure/blob/main/notebooks/exploration.ipynb) | [run](https://colab.research.google.com/github/BobKerns/epi-structure/blob/main/notebooks/exploration.ipynb) | [run](https://mybinder.org/v2/gh/BobKerns/epi-structure/main?labpath=notebooks/exploration.ipynb) |

### Applied Examples

**contact_matrix_exploration.ipynb**: Contact matrix design sandbox
- Explore matrix construction patterns, defaults, and sparse bridge overrides

**exploration.ipynb**: Free-form model exploration
- Scratchpad notebook for ad hoc checks, plots, and parameter experiments

## Contact Matrix UX Pattern

For quick exploratory work, prefer building contact matrices from:

- population names in engine order
- within-group defaults/presets (diagonal entries)
- sparse bridge overrides only where needed

This keeps configuration readable while still producing a full square matrix for the engine.

## Scenario Builder API Pattern

For structured scenarios, prefer the scenario-centered builder API over ad hoc matrix mutation:

- start from a seed or base scenario: `Scenario.builder(base=...)`
- register/adjust populations as needed
- set transitions explicitly with `.transition(recipient, source, value)`
- read or tweak one transition with indexing syntax:
   - lookup: `scenario["general", "cluster"]`
   - assignment: `scenario["general", "cluster"] = new_value`

This keeps scenario definitions and matrix edits explicit, reproducible, and easier to compare across notebook variants.

## Usage

To run the notebooks, ensure you have the required environment set up:

1. Activate the virtual environment:
   ```bash
   source .venv/bin/activate
   ```

2. Start Jupyter Lab:
   ```bash
   jupyter lab
   ```

3. Open the desired notebook from the `notebooks/` directory.

## Notes

- These notebooks are for exploration and demonstration purposes.
- Ensure the `epi_structure` package is installed in your environment before running the notebooks.
