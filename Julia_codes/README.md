# Julia Codes — PhD Research Portfolio

This folder showcases Julia tools I developed during my PhD to study **electronic structure and quantum phenomena in condensed matter systems**.

The repository focuses on reproducible post-processing workflows for:
- **Wannier90 tight-binding models** (`*_hr.dat`, `*_eig`)  
- **VASP outputs** (`PROCAR`, `EIGENVAL`, HDF5-derived datasets)  
- **Semiconductor quantum-dot physics** (electronic coupling and Fano line-shape analyses)

If you are reviewing this portfolio for hiring or collaboration, this section is designed to help you quickly understand what each code does and where it fits in a research workflow.

---

## What is in this folder?

### Main scripts

| File | Research purpose | Typical inputs | Core output |
|---|---|---|---|
| `Bandstructure.jl` | Tight-binding band structures from Wannier Hamiltonians | `seedname_hr.dat`, k-path definitions | Band energies and publication-style plots |
| `DOS.jl` | DOS/PDOS workflows for collinear and non-collinear cases | HDF5 energies and VASP-derived data | DOS/PDOS curves |
| `FS.jl` | 2D/3D Fermi-surface visualization | HDF5 band datasets | Fermi-surface contour maps |
| `TB-Greens.jl` | Surface Green's function via principal layers | Layer-resolved Hamiltonians in HDF5/CSV | Surface spectral maps |
| `VASP_DOS.jl` | PROCAR parsing and DOS/PDOS projections | `PROCAR` files | Total/projection-resolved DOS |
| `VASP_Fermi_surface.jl` | VASP-based Fermi-surface extraction | `EIGENVAL` and converted tabular data | Constant-energy contour plots |
| `Electronic_couplings.jl` | Quantum-dot coupling and wavefunction trends | Quantum-dot size/distance parameter grids | Coupling-vs-size/distance analyses |
| `Fano_IR_1S.jl` | Fano resonance analysis in pump-probe context | Experimental/simulated line-shape arrays | Optimized Fano parameters and trend plots |

### Notebooks

- `TB_Bandstructure.ipynb` — interactive band-structure exploration.  
- `TB_DOS.ipynb` — DOS post-processing from tight-binding data.  
- `TB_FS.ipynb` — interactive Fermi-surface calculations.  
- `VASP_DOS.ipynb` — DOS/PDOS exploration from VASP outputs.

---

## Skills demonstrated

- **Scientific programming in Julia** for medium-to-large post-processing pipelines.
- **Computational condensed matter physics** (band structure, DOS/PDOS, Fermi surfaces).
- **Electronic-structure interoperability** across Wannier90 and VASP ecosystems.
- **Model-driven quantum-dot analysis** (size- and distance-dependent couplings).
- **Data transformation + visualization** using HDF5/CSV and publication-ready plotting.

---

## Quick start

> These scripts were created for research workflows and often include machine-specific paths.  
> For reuse, update file paths and data locations at the top of each script.

1. Install Julia and core packages used in this folder (e.g., `Plots`, `HDF5`, `CSV`, `DataFrames`, `LinearAlgebra`, `DelimitedFiles`, `Distributions`).
2. Open the target `.jl` file and adjust path constants (for example `FILE_DATA_*` or absolute `raw"..."` paths).
3. Run the specific analysis function for your case (see [`FUNCTION_INDEX.md`](./FUNCTION_INDEX.md)).

---

## Suggested reading order (for recruiters)

1. **`Bandstructure.jl` + `TB_Bandstructure.ipynb`**: core tight-binding workflow.  
2. **`VASP_DOS.jl` / `DOS.jl`**: practical parsing + DOS/PDOS analysis depth.  
3. **`FS.jl` + `TB-Greens.jl`**: advanced reciprocal-space and Green-function methods.  
4. **`Electronic_couplings.jl` + `Fano_IR_1S.jl`**: quantum-dot and spectroscopy-focused modeling.

---

## External context

For broader project context and group-level codebases, see the [RibeiroGroup GitHub organization](https://github.com/RibeiroGroup).
