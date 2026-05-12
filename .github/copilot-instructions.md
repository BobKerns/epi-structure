# Repository Instructions

## Implementation Conventions

- Use `uv` for Python package and environment management.
- Prefer Python scripts over shell scripts for project automation and developer workflows.
- Use `click` for command-line parsing in standalone CLI tools.
- For standalone Python scripts, use PEP 723 inline script metadata headers so they can be run directly with `uv run`.
- Keep runtime library code, developer-maintenance tooling, and temporary one-off scripts clearly separated.

## Directory Intent

- `src/epi_structure/`: runtime library code imported by notebooks, applications, or other Python modules.
- `scripts/`: maintained developer utilities and repeatable project automation.
- `tmp/` or another clearly temporary location: disposable one-off experiments created during development. Do not promote temporary scripts into `scripts/` without cleanup.

## Separation Rules

- Do not place ad hoc development scripts inside `src/`.
- Do not put runtime package logic into notebooks or shell entrypoints.
- If a script becomes part of normal project workflow, keep it in `scripts/` and make it clean, documented, and repeatable.
- If logic in a script is needed at runtime, move that logic into `src/epi_structure/` and keep the script as a thin CLI wrapper.

## Package Management

- Prefer `uv sync`, `uv add`, and `uv run` over `pip install`, `python -m venv`, or ad hoc virtualenv commands unless there is a specific reason not to.
- Keep dependency declarations in project configuration, not in scattered setup notes.

## Spellcheck Vocabulary

- Keep shared project vocabulary in the committed `cspell.json` plus `cspell-project-words.txt` files.
- Do not store shared spelling words in `.vscode/settings.json`; reserve that file for truly local, user-specific workspace preferences.

## CLI Style

- Use descriptive command names and explicit options.
- Favor small Python entrypoints that call library functions rather than embedding core logic directly in CLI handlers.
- Keep command output concise and useful for iterative development.

## Temporary Agent Artifacts

- Any temporary files created during implementation should be easy to identify and easy to delete.
- Remove throwaway artifacts before finishing work unless the user explicitly wants them preserved.
