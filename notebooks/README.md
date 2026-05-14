# Notebooks

This directory contains Jupyter notebooks for exploring and demonstrating the epidemic population structure model.

## Current Notebooks

- `exploration.ipynb`: A sample notebook for driving the simulation engine with example data.
- `contact_matrix_exploration.ipynb`: A focused two-population example using `StructuredEpidemicModel` and a lightweight named matrix builder.
- `05_interventions.ipynb`: Illustrates how familiar policy framings map to matrix-scaling interventions, includes a group/bridge link map visualization, and compares scenario outcomes.

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
