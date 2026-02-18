# Julia Codes (PhD Portfolio)

This folder contains Julia scripts used during my PhD for **electronic-structure post-processing**, **tight-binding modeling**, and **atomistic structure manipulation**.
The code demonstrates a practical research workflow that connects first-principles outputs (e.g., VASP data) with fast model analysis and visualization.

---

## What this folder demonstrates

- Building and analyzing **tight-binding Hamiltonians**.
- Computing electronic observables such as **band structure**, **DOS/PDOS**, and **Fermi surfaces**.
- Post-processing and interoperability between **VASP outputs** and Julia analysis scripts.
- Utility workflows for **structure editing** (supercells, coordinate scaling, plane editing, hydrogen addition).

---

## Repository Guide

### Tight-binding and model analysis
- `Bandstructure.jl` – band-structure workflows.
- `Bandstructure_TB.jl` – tight-binding-specific band calculations.
- `TB_PP.jl` – tight-binding post-processing utilities.
- `TB_VASP.jl` – bridge between tight-binding routines and VASP-related data.
- `TB_Fermi_Surface.jl` – tight-binding Fermi-surface evaluation.
- `TB-Greens.jl` – Green's-function-based tight-binding analysis.

### Density of states and projections
- `DOS.jl`, `DOS_m.jl`, `DOS_s.jl` – DOS calculation variants for different analysis contexts.
- `PDOS.jl` – projected density of states workflows.

### VASP-related analysis
- `VASP_Fermi_surface.jl` – Fermi-surface extraction/visualization from VASP-type data.

### Structure manipulation and geometry utilities
- `change_supercell.jl` – supercell transformations.
- `scale_atoms.jl` – atomic coordinate/lattice scaling.
- `chang_atoms.jl` – atomic-position editing utility.
- `fix_plane.jl`, `plane_eraser.jl` – plane-level geometry cleanup/editing.
- `Adding_hydrogens.jl` – hydrogen passivation/addition routine.
- `crystal_plot.jl` – crystal structure plotting helper.

### Miscellaneous
- `Old_codes.jl` – archived/legacy routines.
- `POSCAR` – example structure input used by some scripts.

---

## How to run

1. Open a terminal in `Julia_codes/`.
2. Run a script with Julia, for example:
   ```bash
   julia Bandstructure_TB.jl
   ```
3. For scripts that depend on external files (e.g., `POSCAR` or VASP outputs), keep required inputs in the expected relative paths.

---

## Suggested reading order (for recruiters/collaborators)

1. `TB_VASP.jl` and `Bandstructure_TB.jl` (core model + workflow integration).
2. `DOS.jl` / `PDOS.jl` (electronic observables).
3. `VASP_Fermi_surface.jl` and `TB_Fermi_Surface.jl` (advanced reciprocal-space analysis).
4. Structure-editing utilities (`change_supercell.jl`, `scale_atoms.jl`, `Adding_hydrogens.jl`) to see practical data-engineering support for simulations.

---

## Why this is useful for hiring review

This folder highlights transferable skills for computational physics/materials and quantitative R&D roles:

- converting physics ideas into robust, scriptable computation,
- integrating ab initio data with reduced models,
- creating reproducible analysis pipelines with clear utility tooling.
