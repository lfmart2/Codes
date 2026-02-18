# MATLAB Codes for Condensed Matter & Molecular Dynamics

This folder collects MATLAB scripts developed during my PhD work on **topological materials**, **tight-binding models**, and **electron–vibration effects**.
The goal of this section is to make the research workflow easy to browse, reproduce, and extend.

---

## Research Focus

- **Topological band theory** (Wilson loops, SSH model physics).
- **Electron occupancy and charge transfer** in adsorbate–polyacetylene systems.
- **Electronic friction** from vibronic coupling (non-adiabatic correction to Born–Oppenheimer dynamics).
- **Data generation for topological defects** (soliton-like structures in staggered chains).

---

## Repository Map (MATLAB)

| File | Type | What it does | Typical output |
|---|---|---|---|
| `Wilson_Loop.m` | Script | Computes Wilson-loop eigenphases for a 4×4 Bloch Hamiltonian over a 2D Brillouin zone and plots Wannier-sector flow. | Wilson-loop phase spectra / Wannier bands |
| `ElectronOccupancy.m` | Script | Studies electron donation/occupancy in SSH-like chains with adsorbates at edge and bulk positions; exports CSV data for plotting/analysis. | Occupancy curves and `.csv` datasets |
| `SSH_bulk.mlx` | Live Script | Bulk SSH model exploration (topological/trivial regimes and spectra). | Interactive figures and parameter scans |
| `SSH_edge.mlx` | Live Script | Edge-state analysis in finite SSH chains. | Edge-localized state plots |
| `El_Friction.mlx` | Live Script | Electronic friction tensor workflow from vibronic coupling setup. | Friction-related observables/figures |
| `El_Occ_Data_Generator.mlx` | Live Script | Generates datasets for topological defects (soliton center/width control). | Defect-resolved datasets |

---

## How to Run

1. Open MATLAB in this folder (or add this folder to the MATLAB path).
2. Start with one of the two main scripts:
   - `Wilson_Loop.m`
   - `ElectronOccupancy.m`
3. For step-by-step visual exploration, open the `.mlx` files in the Live Editor.

> Notes:
> - Some scripts assume helper functions (e.g., SSH chain builders) are available in your MATLAB path.
> - Numerical resolution parameters (`N`, k-point grids, coupling constants) can be increased for publication-quality convergence.

---

## Suggested Reading Order (for recruiters/collaborators)

If you want a fast overview of technical depth:

1. **`Wilson_Loop.m`** → topological invariant workflow and Brillouin-zone numerics.
2. **`ElectronOccupancy.m`** → model building, parameter sweeps, physical interpretation, and data export.
3. **`.mlx` live scripts** → interactive derivations, diagnostics, and presentation-ready figures.

---

## Why this folder matters

These scripts reflect practical PhD-level skills in:

- Translating condensed-matter theory into robust numerical pipelines.
- Designing reproducible parameter scans and post-processing workflows.
- Connecting physical observables (occupancy, Wilson phases, friction) to interpretable figures and datasets.
