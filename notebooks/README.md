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

**02_seir.ipynb** (upcoming): SEIR model adds latent period
- S → E → I → R progression
- How latent period affects timing and peak while R₀ remains constant
- Comparison with SIR baseline

**03_seirs.ipynb** (upcoming): SEIRS adds waning immunity
- Full loop with loss of immunity
- Endemic equilibrium and cycle behavior
- How R₀ relates to endemic prevalence

**04_two_population.ipynb** (upcoming): Multi-population dynamics
- Two populations (general + cluster) with contact matrix
- Within-group and bridge transmission
- How coupling affects each population differently
- Effective R₀ is population-dependent

### Applied Examples

- **05_interventions.ipynb**: Policy language to matrix-scaling transformations
  - Ventilation, tracing, distancing as transmission scalar interventions
  - Group/bridge link visualization with subset-selection + transformation pattern
  - Scenario comparison (baseline, cluster ventilation, bridge tracing, targeted combo)

## Contact Matrix UX Pattern

For quick exploratory work, prefer building contact matrices from:

- population names in engine order
- within-group defaults/presets (diagonal entries)
- sparse bridge overrides only where needed

This keeps configuration readable while still producing a full square matrix for the engine.

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
