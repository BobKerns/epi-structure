# Epidemic Population Structure: Interactive Demonstration
## Design Document

*Captures design decisions and context from initial brainstorming session, May 2026.*

---

## Project Goal

An interactive demonstration showing how population structure — overlapping subpopulations, superspreader clusters, and bridge individuals — drives epidemic dynamics in ways that homogeneous models miss. The project has three purposes:

1. **Educational** — build intuition from first principles, accessible to people with no prior SIR knowledge
2. **Explanatory** — show why extrapolation between populations is dangerous
3. **Advocacy** — demonstrate concretely how targeted interventions work, in support of sound public health policy

---

## Motivation and Context

Arising from a Twitter thread (May 7, 2026) explaining why population structure must be taken into account when assessing epidemic risk, and why superspreader events connected to other superspreader events change the calculus fundamentally. The thread used a diagram (generated with Claude) showing isolated vs. bridge-connected clusters.

The key policy insight — that targeted intervention in clusters and bridges can collapse an epidemic without acting on the general population at all — is underappreciated and has direct relevance to current public health debates.

---

## The Lynchpin Insight

The conceptual heart of the demonstration:

> A general population with R₀ < 1 **cannot sustain an epidemic on its own**, regardless of how many times it is seeded. This is qualitatively different from susceptible depletion. Add a supercritical cluster (R₀ > 1) with a bridge into the general population, and continuous reseeding keeps the epidemic alive. Cut the bridge — or reduce cluster R₀ below 1 — and the general population clears, **without changing its own R₀ at all**.

This is the moment the demonstration must make visceral, not just intellectual.

---

## Pedagogical Ladder

Complexity accretes one concept at a time. Each stage adds exactly one new idea, and the visualization reflects only what has been introduced:

1. **Single population, SIR** — the basic epidemic wave. R₀ as sole parameter. Outbreak dies even without intervention (susceptible depletion). Establish R₀ = 1 as the critical threshold.
2. **Add E compartment → SEIR** — latent period smooths and delays the curve. Why early detection is hard.
3. **Add waning immunity → SEIRS** — endemic equilibria, oscillations. Why Covid didn't just "burn out."
4. **Two populations, no bridge** — general population (R₀ < 1) plus one cluster (R₀ > 1). Cluster burns locally, general population unaffected. Feel the certainty of subcritical extinction.
5. **Add the bridge** — continuous reseeding. The lynchpin moment.
6. **Interventions** — reduce cluster R₀ (ventilation, masking), cut bridges (contact tracing, isolation), targeted vaccination. Each is a slider or toggle. Viewer sees that general population R₀ need not change.
7. **Multiple clusters, varying bridge rates** — full model.

UI design implication: controls for later stages are hidden or grayed out until earned by the narrative progression.

---

## Emotional/Narrative Arc

- Seed the subcritical general population repeatedly — it always dies out. Feel the certainty.
- Add an isolated supercritical cluster — burns locally, general population unaffected.
- Add the bridge — general population now gets continuously reseeded.
- Cut the bridge — general population clears even without changing its own R₀.
- **Policy punchline**: we had targeted options. They required infrastructure and political will. We should build and protect that infrastructure now.

---

## Model

**SEIRS with contact matrix** — small number of subpopulations (general population + 2–3 cluster types), with a matrix of between-group contact rates encoding bridge structure.

**Numerical integration**: simple iterative RK4 (Runge-Kutta 4th order). Chosen over scipy.integrate.solve_ivp because:
- Integration loop is readable as the whiteboard equations
- Can pause, inspect state at any timestep, animate step-by-step
- Animation *is* the integration loop, not a replay of precomputed results
- Adequate for epidemiologically realistic parameter ranges

Stepsize exposed as an inspectable parameter.

---

## Implementation Architecture

### Python package as canonical implementation

Simulation kernel lives in a proper Python package (working name: `epi_structure`). Single source of truth. Clean API. Version controlled.

**Minimal initial structure:**
```
epi_structure/
    src/
        epi_structure/
            __init__.py
            model.py        # SEIRS kernel, contact matrix, RK4 integration
            parameters.py   # Parameter dataclasses
            interventions.py # Intervention logic
    notebooks/
        01_sir_single.ipynb
        02_seir.ipynb
        03_seirs.ipynb
        04_two_populations.ipynb
        05_interventions.ipynb
    DESIGN.md
    README.md
    pyproject.toml
```

No tests, CI, or full documentation initially — just enough structure that notebooks import from the package from day one.

### Jupyter notebooks

Import from the package. Focus entirely on narrative and visualization. No simulation code in the notebook cells themselves. ipywidgets for sliders. Development and family-accessible exposition environment.

### Public demo (later)

**Pyodide** (Python running in browser via WebAssembly) is the preferred bridge — the actual Python package potentially runs client-side with no server required. Observable HQ as the public-facing platform (account already exists). Investigate Pyodide viability before committing to a JavaScript port.

Shared markdown narrative text, dynamically imported in Observable, keeps exposition in sync between notebook and public demo as the project iterates.

---

## Deployment Targets

| Artifact | Audience | Format |
|---|---|---|
| Jupyter notebooks | Family, technically-minded readers, those who want to see the code | Binder or Google Colab link |
| Observable notebook | General public | Shareable URL, no setup |
| Journal paper (future) | Academic credibility, peer review | PLOS Computational Biology or Epidemics |

---

## Phase 1 Definition (Current Build Target)

Phase 1 is a notebook-first, research-grade MVP focused on making the lynchpin insight clear and testable.
It prioritizes correctness, inspectability, and narrative flow over production polish.

### Primary Objective

Deliver a working interactive demonstration that proves, in simulation, all of the following:

1. A subcritical general population (R0 < 1) clears repeatedly seeded outbreaks.
2. A supercritical cluster (R0 > 1) can sustain local transmission.
3. Bridge-driven reseeding from cluster to general population can maintain ongoing general-population incidence.
4. Targeted interventions (cut bridge or reduce cluster R0) clear the general population without changing general-population R0.

### In-Scope Artifacts

1. Python package with SEIRS + contact-matrix kernel and RK4 integration.
2. Notebook sequence through pedagogical stages 1 to 6.
3. Interactive controls in notebooks (ipywidgets) for parameters and interventions.
4. Reproducible baseline parameter sets for each stage.
5. Static plots and time-series outputs suitable for screenshots/shareable discussion.
6. Minimal README with setup and execution instructions.

### Explicitly Out of Scope (Phase 1)

1. Public Observable deployment.
2. Pyodide browser execution.
3. Multi-cluster full model from stage 7.
4. Production UI/branding work.
5. Formal manuscript drafting.
6. Full CI/release pipeline.

### Functional Requirements

1. Deterministic integration for fixed parameters and seed conditions.
2. Per-group compartments S, E, I, R and optional waning immunity (R -> S).
3. Configurable contact matrix with directional between-group transmission.
4. Toggleable intervention primitives:
    - Scale within-cluster transmission.
    - Scale bridge transmission terms.
    - Apply at configurable simulation times.
5. Outputs at each timestep for each compartment and group.
6. Notebook controls that do not require code editing to run core scenarios.

### Narrative Requirements

1. Each notebook introduces one new concept only.
2. The same baseline parameters are reused where possible to reduce cognitive churn.
3. Every stage includes one sentence answering: "What did this stage add that the prior stage could not explain?"
4. Stage 6 includes a dedicated "policy interpretation" cell tied directly to intervention controls.

### Milestones and Exit Criteria

1. M1: Package skeleton and imports working in notebook.
    Exit: notebook can import package and run one SIR trajectory end-to-end.
2. M2: Single-population stages 1 to 3 complete.
    Exit: SIR, SEIR, and SEIRS figures generated from package-only simulation code.
3. M3: Two-population no-bridge and with-bridge scenarios complete.
    Exit: side-by-side plots show extinction vs sustained reseeding behavior.
4. M4: Intervention controls wired for stage 6.
    Exit: at least two interventions demonstrably clear general-population incidence while general-population R0 is unchanged.
5. M5: Documentation pass.
    Exit: README run instructions + parameter table + known limitations section are present.

### Validation Scenarios (Must Pass)

1. Subcritical-alone check:
    General population only, R0 < 1, repeated seeding pulses -> infections decay after each pulse.
2. Isolated cluster check:
    Cluster R0 > 1, no bridge -> general population remains near baseline infection level.
3. Bridge persistence check:
    Cluster R0 > 1 with bridge -> general population shows continued incidence.
4. Bridge cut check:
    Disable bridge mid-run -> general-population incidence decays without changing its own R0.
5. Cluster control check:
    Reduce cluster R0 below 1 mid-run -> reseeding collapses and general population clears.

### Parameter Strategy for Phase 1

1. Start with normalized populations (for example, each group size = 1.0) to make structure effects visually clear.
2. Use a shared default latent and infectious period unless a stage explicitly demonstrates timing effects.
3. Keep intervention strengths coarse (for example, 0%, 25%, 50%, 75%, 100%) for interpretability.
4. Document all baseline values in a single table reused across notebooks.

### Risks and Mitigations

1. Risk: parameter choices accidentally obscure the core effect.
    Mitigation: lock one canonical baseline set and test sensitivity with small perturbations.
2. Risk: notebook code drifts from package logic.
    Mitigation: no simulation equations in notebooks; package imports only.
3. Risk: over-complicated controls reduce clarity.
    Mitigation: hide advanced controls behind an "advanced" toggle.
4. Risk: policy claim appears stronger than model scope supports.
    Mitigation: include a "model limits" callout in every stage from 4 onward.

### Phase 1 Deliverable Checklist

1. Package modules: model, parameters, interventions.
2. Notebooks 01 through 05 runnable from clean environment.
3. Parameter table and scenario definitions documented.
4. At least one exportable figure per stage.
5. Written limitations and non-claims section.

---

## Collaborators

**Erika Kerns** (coauthor) — MS Disease Ecology, UCLA (completing 2026). Thesis: hierarchical Bayesian modeling of leptospirosis population dynamics in California sea lions. Background in necropsy, field sampling, TMMC and Channel Islands data. Entering Colorado State veterinary school. Provides domain expertise and academic credibility. Whiteboard origin story (2016) is a natural narrative hook.

---

## Key Design Principles

- **Inspectability over black-box results** — the process of the simulation is as important as the outcome
- **Accretive complexity** — never introduce two concepts simultaneously
- **Implement once** — no backend server; Python package is canonical; Pyodide bridges to public demo
- **Library from the start** — minimal package structure before notebooks, not extracted afterward
- **Policy advocacy is the end goal** — the demonstration exists to make a concrete, actionable argument, not just to educate

---

## Open Questions

- Pyodide feasibility for the Observable deployment (investigate early)
- Exact cluster types to include in initial demo (party/transient, dorm/persistent, general population confirmed; cruise ship optional)
- Whether shared markdown narrative text approach works cleanly in Observable
- Journal paper timing and division of authorship responsibilities

---

*This document should be updated as design decisions evolve. It doubles as context for AI-assisted development sessions.*
