# Diagnostics Hooks

This document outlines potential diagnostic hooks for the genetic algorithm (GA) and simulation routines.

## GA Progress
- Use a custom output function in `viscous.m` to record each generation's best score and constraints.
- `ga.opt.save_log` can be extended to include intermediate population statistics for live monitoring.

## Jacobian Sparsity
- Capture the Jacobian structure from `mck_with_damper_adv` inside `simulate.m` and store it for visualization.
- MATLAB/Octave's `odeset` with `JPattern` can expose sparsity patterns when solving the ODEs.

## Constraint Details
- Within the GA loop in `viscous.m`, log the active constraints and penalties for each candidate.
- `simulate.m` already returns `resp.metrics`; expanding this to include constraint violation details can aid debugging.

## Integration Points
- **simulate.m**: insert hooks around the call to `mck_with_damper_adv` to log Jacobian patterns and constraint evaluations.
- **viscous.m**: add GA output functions after configuration of `ga.opt.*` to capture progress and constraint summaries per generation.

These hooks can be enabled conditionally to avoid performance penalties during routine runs.
