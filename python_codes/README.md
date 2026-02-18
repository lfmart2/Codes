# Python Codes — PhD Research Portfolio

This folder contains Python workflows developed during my PhD for **electronic-structure post-processing and local density of states (LDOS) analysis**.

The main focus is high-throughput simulation tooling around:
- **Wannier90 real-space Hamiltonians** (`wannier90_hr.dat`)  
- **Kwant-based tight-binding supercell construction**  
- **DOS/LDOS and hydrogen-coupling analysis** in SOC and non-SOC regimes  
- **Optional GPU acceleration (CuPy)** for selected LDOS kernels

---

## Folder map

### Main scripts

| File | Purpose | Typical inputs | Output |
|---|---|---|---|
| `DOS_kwant_nSOC_supercell.py` | Non-SOC DOS pipeline from Wannier90 + Kwant supercell | `wannier90_hr.dat`, hydrogen coupling files | DOS arrays / serialized results |
| `DOS_kwant_SOC_supercell.py` | SOC DOS pipeline with spinful Hamiltonian handling | `wannier90_hr.dat`, hydrogen coupling files | DOS arrays / serialized results |
| `LDOS_kwnat_nSOC_supercell.py` | CPU LDOS workflow (legacy filename spelling retained) | Wannier + coupling inputs | LDOS and projected spectra |
| `LDOS_kwnat_SOC_supercell.py` | CPU LDOS workflow with SOC (legacy filename spelling retained) | Wannier + coupling inputs | LDOS and projected spectra |
| `LDOS_kwant_nSOC_supercell_gpu.py` | nSOC LDOS workflow with optional GPU arrays | Same as nSOC CPU script | Faster LDOS kernels with `--use-gpu` |
| `LDOS_kwant_SOC_supercell_gpu.py` | SOC LDOS workflow with optional GPU arrays | Same as SOC CPU script | Faster LDOS kernels with `--use-gpu` |
| `LDOS_gpu_test.py` | Experimental GPU-integration prototype utilities | Coupling/model data | Benchmark/prototype LDOS outputs |

### Notebooks

- `Electron_Occupancy.ipynb` — occupancy-oriented exploration notebook.  
- `Electron_Occupancy_solitons.ipynb` — occupancy/soliton analysis notebook.

---

## Skills demonstrated

- **Scientific Python engineering** (`numpy`, `h5py`, `pandas`) in research workflows.
- **Quantum transport / tight-binding modeling** with `kwant`.
- **Model parsing and data-pipeline design** from Wannier90-derived files.
- **Performance-minded implementation**, including optional **GPU acceleration with CuPy**.
- **SOC vs non-SOC workflow handling** for realistic materials simulations.

---

## Quick start

> These scripts are research-grade and may include machine-specific paths/constants.  
> Update paths and runtime parameters before first execution.

1. Create a Python environment with required packages (typical: `numpy`, `scipy`, `pandas`, `h5py`, `matplotlib`, `kwant`; optionally `cupy`).
2. Inspect the target script and update input/output paths and supercell settings.
3. Run a CPU baseline first, then enable GPU mode where available.

GPU example commands:

```bash
python3.10 LDOS_kwant_nSOC_supercell_gpu.py \
  <input_tb_up> <output_up> \
  <input_tb_dn> <output_dn> \
  --use-gpu
```

```bash
python3.10 LDOS_kwant_SOC_supercell_gpu.py \
  <input_tb_soc> <output_soc> \
  --use-gpu
```

For a function-level navigation guide, see [`FUNCTION_INDEX.md`](./FUNCTION_INDEX.md).

---

## Suggested reading order (for recruiters)

1. `DOS_kwant_nSOC_supercell.py` and `DOS_kwant_SOC_supercell.py` (core architecture).  
2. `LDOS_kwant_*_gpu.py` (performance-aware implementation + GPU path).  
3. `LDOS_gpu_test.py` (prototype/experimental acceleration work).  
4. Notebooks for exploratory, analysis-facing outputs.
