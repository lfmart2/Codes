# Function Index (Julia Codes)

This quick index maps the main entry points in each Julia script to their scientific purpose.

## `Bandstructure.jl`

- `load_hr(path)` — parses Wannier90 `*_hr.dat` Hamiltonians.
- `bloch_hamiltonian(hr, k)` — builds the Bloch Hamiltonian at a k-point.
- `tb_proj_bands_*` family — computes tight-binding projected bands for SOC/non-SOC variants.
- `bandstructure_*` family — generates band plots for multiple spin/orbital scenarios.
- `TBvsVASP_*` family — compares tight-binding and VASP bands.

## `DOS.jl`

- `DOS_nSOC(...)` / `DOS_SOC(...)` / `DOS_SOC_H(...)` — DOS broadening workflows.
- `DOS_PROCAR_*` family — DOS from VASP-derived `PROCAR` data.
- `PDOS_PROCAR_*` family — projection-resolved DOS (including H-passivated variants).

## `FS.jl`

- `plot_FS_collinear(...)`, `plot_FS_ncollinear(...)` and related functions — Fermi-surface contour generation for different spin treatments.
- `test_FS(...)`, `test_FS_COHP(...)` — analysis/testing helpers for FS generation.

## `TB-Greens.jl`

- `load_planes(...)` / `load_planes_int(...)` — layer data loaders.
- `H_k(...)`, `T_k(...)` — intra-/inter-layer k-space Hamiltonian blocks.
- `bulk_Greens(...)`, `surf_Greens(...)` — Green's function solvers.
- `plot_FS(...)` — Fermi-surface-like spectral visualization from Green's functions.

## `VASP_DOS.jl`

- `load_PROCAR_*` family — parsers for DOS/PDOS/magnetization from `PROCAR`.
- `DOS_PROCAR_*` and `PDOS_PROCAR_*` families — plotting/analysis for SOC and non-SOC runs.
- `FD_dist(...)` — Fermi-Dirac helper.

## `VASP_Fermi_surface.jl`

- `filter_file(...)` — cleans/parses tabularized VASP band information.
- `Fermi_surface(...)` — computes and plots constant-energy surfaces from VASP outputs.

## `Electronic_couplings.jl`

- `cond_bnd_edge(...)` — conduction-band edge model helper.
- `t_size_*` / `t_dis_*` families — coupling vs size/distance models.
- `ψ_size_plot_*` / `ψ_dis_plot_*` families — wavefunction/coupling visual summaries.

## `Fano_IR_1S.jl`

- `qparam(...)` — Fano q-parameter model helper.
- `sse_1S(...)` — objective function.
- `optimize_params_constrained_1S(...)` — constrained fitting routine.
- `qq_ds_Fano_1S()` / `qq_sz_Fano_1S()` — fitted trend analyses.

---

If you are new to this folder, start with `Bandstructure.jl`, `VASP_DOS.jl`, and the associated notebooks for a fast overview of both workflow maturity and physics coverage.
