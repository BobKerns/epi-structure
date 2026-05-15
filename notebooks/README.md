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

**04_two_population.ipynb**: Multi-population dynamics
- Two populations (general + cluster) with contact matrix
- Within-group and bridge transmission
- How coupling affects each population differently
- Effective R₀ is population-dependent



**05_superspreader_party_crawl.ipynb**: Sparse coupling and re-seeding dynamics
- General-only subcritical regime ($R_0 < 1$) contrasted with coupled cluster dynamics
- Superspreader cluster surges and spillback into the general population

**06_rt_wells_riley_decomposition.ipynb**: Mechanistic $R_t$ decomposition
- Decomposes $R_t$ into contact rate, infectious duration, susceptibility, and per-contact infection probability
- Uses a Wells-Riley style dose model for infection probability terms

**07_interventions.ipynb**: Policy language to matrix-scaling transformations
   - Ventilation, tracing, distancing as transmission scalar interventions
   - Group/bridge link visualization with subset-selection + transformation pattern
   - Scenario comparison (baseline, cluster ventilation, bridge tracing, targeted combo)

**08_timed_interventions.ipynb**: Time-varying interventions
- Component-level policy schedules for $c_t$ and Wells-Riley terms
- Time-bounded matrix interventions in the structured model

### Applied Examples

-- **07_interventions.ipynb**: Policy language to matrix-scaling transformations
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
