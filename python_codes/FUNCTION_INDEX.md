# Function Index (Python Codes)

This document summarizes the main callable entry points in the Python scripts under `python_codes/`.

## `DOS_kwant_nSOC_supercell.py`

- `HRData` — structured container for parsed Wannier90 real-space Hamiltonian data.
- `load_hr(path)` — parser for `wannier90_hr.dat` into `HRData`.
- `_build_hydrogen_coupling_map()` — precomputes label-to-index mapping for hydrogen interactions.
- `coupling_hydrogen_slab_mapping(label)` — returns coupling indices for a given interaction label.
- `_load_hydrogen_block_keys(file_path, r_vec)` — reads available interaction keys from HDF5 block.
- `coupling_hydrogen_slab(...)` — constructs hydrogen-coupling matrix blocks.
- `iter_hoppings_with_hydrogen(...)` — iterates slab hoppings plus hydrogen couplings.
- `build_supercell_with_hydrogen(...)` — assembles Kwant supercell model.
- `report_system_size(fsys)` — helper to print/report final system size.

## `DOS_kwant_SOC_supercell.py`

- Mirrors the non-SOC architecture for SOC-enabled models.
- Same key components: `HRData`, `load_hr`, coupling map helpers, hopping iterators, supercell builder, and system-size reporting.

## `LDOS_kwnat_nSOC_supercell.py` (CPU)

- `HRData`, `load_hr`, and hydrogen-coupling helper functions for CPU LDOS runs.
- `compute_hydrogen_pdos_kpm(...)` — computes projected DOS/LDOS with KPM-based flow.
- `quick_test_hydrogen_pdos(...)` — convenience test driver for quick validation.

## `LDOS_kwnat_SOC_supercell.py` (CPU)

- SOC counterpart to the non-SOC CPU LDOS workflow.
- Includes parsing, coupling assembly, supercell building, and KPM PDOS evaluation helpers.

## `LDOS_kwant_nSOC_supercell_gpu.py`

- `_get_array_module(use_gpu)` — switches NumPy/CuPy backend.
- `_to_numpy(array)` — normalizes arrays back to CPU memory when needed.
- `compute_hydrogen_pdos_kpm(..., use_gpu=False)` — GPU-aware PDOS/LDOS kernel path.
- `quick_test_hydrogen_pdos(...)` — fast validation entry point.
- Includes the same HR parser and coupling-construction functions used in CPU workflows.

## `LDOS_kwant_SOC_supercell_gpu.py`

- SOC GPU variant with the same backend-switch pattern and helper architecture.
- Main workflow is centered around GPU-aware coupling assembly and PDOS/LDOS evaluation.

## `LDOS_gpu_test.py`

- `GPUEnhancedKwantKPM` — experimental class for integrated GPU workflows.
- `compute_hydrogen_pdos_kpm_gpu_integrated(...)` — prototype integrated GPU PDOS routine.
- `build_hydrogen_nbp_supercell_optimized(...)` — optimized model construction helper.
- `_preload_coupling_data_gpu(...)` and `_compute_hydrogen_coupling_gpu(...)` — low-level GPU utility methods.

---

If you are reviewing this folder quickly, start with:
1) `DOS_kwant_nSOC_supercell.py`, 2) `LDOS_kwant_nSOC_supercell_gpu.py`, and 3) `LDOS_kwant_SOC_supercell_gpu.py`.
