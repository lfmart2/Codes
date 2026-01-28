from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import h5py
import kwant
import numpy as np


@dataclass(slots=True)
class HRData:
    num_wann: int
    nrpts: int
    R: np.ndarray
    weight: np.ndarray
    H_R: np.ndarray


def load_hr(path: str | Path) -> HRData:
    """Load Wannier90 real-space Hamiltonian (wannier90_hr.dat).

    Returns
    -------
    HRData
        Parsed Wannier90 real-space Hamiltonian data. The arrays are:
        - R: integer lattice vectors with shape (3, nrpts)
        - weight: integer weights with shape (nrpts,)
        - H_R: complex hopping blocks with shape (nw, nw, nrpts)
    """
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        handle.readline()  # header
        num_wann = int(handle.readline())
        nrpts = int(handle.readline())

        weight: list[int] = []
        while len(weight) < nrpts:
            line = handle.readline()
            if not line:
                break
            line = line.strip()
            if line:
                weight.extend(int(x) for x in line.split())
        weight_arr = np.array(weight[:nrpts], dtype=int)

        R = np.zeros((3, nrpts), dtype=int)
        H_R = np.zeros((num_wann, num_wann, nrpts), dtype=np.complex128)

        R_index_map: dict[tuple[int, int, int], int] = {}
        idx_counter = 0

        for line in handle:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 7:
                continue

            Rx, Ry, Rz = map(int, parts[0:3])
            i, j = map(int, parts[3:5])
            real, imag = map(float, parts[5:7])

            key = (Rx, Ry, Rz)
            idx = R_index_map.get(key)
            if idx is None:
                if idx_counter >= nrpts:
                    raise ValueError(
                        "Parsed more R points than declared in header; "
                        "check the Wannier90 file."
                    )
                idx = idx_counter
                R_index_map[key] = idx
                R[:, idx] = [Rx, Ry, Rz]
                idx_counter += 1

            H_R[i - 1, j - 1, idx] = real + 1j * imag

        nr_eff = idx_counter
        return HRData(
            num_wann=num_wann,
            nrpts=nr_eff,
            R=R[:, :nr_eff],
            weight=weight_arr[:nr_eff],
            H_R=H_R[:, :, :nr_eff],
        )


def _build_hydrogen_coupling_map() -> dict[str, np.ndarray]:
    mapping: dict[str, np.ndarray] = {}
    layer_specs = {
        "first_d_layer": range(0, 5),
        "second_d_layer": range(5, 10),
        "third_d_layer": range(10, 15),
        "first_sp_layer": range(140, 144),
        "second_sp_layer": range(144, 148),
        "third_sp_layer": range(148, 152),
    }
    for prefix, interaction_range in layer_specs.items():
        interaction_idx = np.fromiter(interaction_range, dtype=int)
        for suffix in ("", "_im"):
            label = f"{prefix}{suffix}"
            mapping[label] = interaction_idx
    return mapping


HYDROGEN_COUPLING_MAP = _build_hydrogen_coupling_map()


def coupling_hydrogen_slab_mapping(label: str) -> np.ndarray | None:
    """Map hydrogen coupling labels to integer indices."""
    return HYDROGEN_COUPLING_MAP.get(label)


@lru_cache(maxsize=32)
def _load_hydrogen_block_keys(
    file_path: str, r_vec: tuple[int, int, int]
) -> list[str]:
    with h5py.File(file_path, "r") as file:
        block_interaction = file[f"[{r_vec[0]},{r_vec[1]},{r_vec[2]}]"]
        return list(block_interaction.keys())


def coupling_hydrogen_slab(
    R_vec: tuple[int, int, int],
    num_slab: int,
    mol_distance: np.ndarray,
    file_path: str | Path,
    x: int,
    y: int,
    Lx: int,
    Ly: int,
    mol_sites: dict[tuple[int, int], int],
) -> np.ndarray:
    """Compute hydrogen coupling block for a given lattice vector."""
    hydrogen_orbitals = 1
    final_block = np.zeros(
        (num_slab + hydrogen_orbitals, num_slab + hydrogen_orbitals), dtype=complex
    )
    rows = np.array([0], dtype=int)
    src_idx = mol_sites.get((x, y))
    dst_idx = mol_sites.get(((x + R_vec[0]) % Lx, (y + R_vec[1]) % Ly))
    has_src = src_idx is not None
    has_dst = dst_idx is not None
    if not has_src and not has_dst:
        return final_block
    distances = (
        np.array([mol_distance[src_idx]], dtype=float) if has_src else np.array([])
    )
    distances_conj = (
        np.array([mol_distance[dst_idx]], dtype=float) if has_dst else np.array([])
    )

    list_Ri = _load_hydrogen_block_keys(str(file_path), R_vec)
    with h5py.File(file_path, "r") as file:
        block_interaction = file[f"[{R_vec[0]},{R_vec[1]},{R_vec[2]}]"]
        for label in list_Ri:
            if label.endswith("_im"):
                continue
            interaction_indices = coupling_hydrogen_slab_mapping(label)
            if interaction_indices is None:
                continue
            data_re = block_interaction[label]
            data_im = block_interaction[f"{label}_im"]
            cols = interaction_indices + hydrogen_orbitals
            if has_src:
                values = (
                    distances[:, None] * data_re[:, 1][None, :]
                    + data_re[:, 0][None, :]
                ) + 1j * (
                    distances[:, None] * data_im[:, 1][None, :]
                    + data_im[:, 0][None, :]
                )
                final_block[np.ix_(rows, cols)] += values
            if has_dst:
                values_conj = (
                    distances_conj[:, None] * data_re[:, 1][None, :]
                    + data_re[:, 0][None, :]
                ) - 1j * (
                    distances_conj[:, None] * data_im[:, 1][None, :]
                    + data_im[:, 0][None, :]
                )
                final_block[np.ix_(cols, rows)] += values_conj.T
    return final_block


def _rvec_is_positive(r_vec: tuple[int, int, int]) -> bool:
    Rx, Ry, Rz = r_vec
    return (Rz > 0) or (Rz == 0 and Ry > 0) or (Rz == 0 and Ry == 0 and Rx > 0)


def iter_hoppings_with_hydrogen(
    hr: HRData,
    mol_distance: np.ndarray,
    coupling_file: str | Path,
    x: int,
    y: int,
    Lx: int,
    Ly: int,
    mol_sites: dict[tuple[int, int], int],
) -> Iterable[tuple[tuple[int, int, int], np.ndarray]]:
    """Yield non-onsite hopping blocks with hydrogen coupling added."""
    hydrogen_orbitals = 1
    block = np.zeros(
        (hr.num_wann + hydrogen_orbitals, hr.num_wann + hydrogen_orbitals),
        dtype=complex,
    )
    for idx in range(hr.nrpts):
        R_vec = tuple(int(v) for v in hr.R[:, idx])
        if R_vec == (0, 0, 0):
            continue
        if not _rvec_is_positive(R_vec):
            continue
        if R_vec[2] != 0:
            continue
        block[:, :] = 0.0
        block[hydrogen_orbitals:, hydrogen_orbitals:] = (
            hr.H_R[:, :, idx] / hr.weight[idx]
        )
        dst_site = ((x + R_vec[0]) % Lx, (y + R_vec[1]) % Ly)
        if (x, y) not in mol_sites and dst_site not in mol_sites:
            yield R_vec, block
            continue
        hydrogen_block = coupling_hydrogen_slab(
            R_vec,
            num_slab=hr.num_wann,
            mol_distance=mol_distance,
            file_path=coupling_file,
            x=x,
            y=y,
            Lx=Lx,
            Ly=Ly,
            mol_sites=mol_sites,
        )
        yield R_vec, block + hydrogen_block


def build_supercell_with_hydrogen(
    hr: HRData,
    Lx: int,
    Ly: int,
    coupling_file: str | Path,
    poly_coeffs_up: list[float],
    seed: int | None = None,
) -> kwant.system.FiniteSystem:
    """Build a 2D supercell including hydrogen coupling data (non-SOC)."""
    num_molecules = (Lx * Ly) // 2
    hydrogen_orbitals = 1
    norb = hr.num_wann + hydrogen_orbitals
    lat = kwant.lattice.square(norbs=norb)
    syst = kwant.Builder()

    rng = np.random.default_rng(seed)
    mol_distance = rng.uniform(-0.6, 0.6, size=hydrogen_orbitals * num_molecules)

    onsite_block = np.zeros((norb, norb), dtype=complex)
    for idx in range(hr.nrpts):
        if np.all(hr.R[:, idx] == 0):
            onsite_block[hydrogen_orbitals:, hydrogen_orbitals:] = hr.H_R[
                :, :, idx
            ].copy()
            break

    mol_sites: dict[tuple[int, int], int] = {}
    mol_idx = 0
    for x in range(Lx):
        for y in range(Ly):
            if (x + y) % 2 != 0:
                continue
            if mol_idx >= mol_distance.size:
                continue
            mol_sites[(x, y)] = mol_idx
            mol_idx += hydrogen_orbitals

    for x in range(Lx):
        for y in range(Ly):
            onsite = onsite_block.copy()
            mol_idx = mol_sites.get((x, y))
            if mol_idx is not None:
                onsite[0, 0] = np.polyval(poly_coeffs_up, mol_distance[mol_idx])
            syst[lat(x, y)] = onsite

    for x in range(Lx):
        for y in range(Ly):
            for (Rx, Ry, Rz), hop_block in iter_hoppings_with_hydrogen(
                hr,
                mol_distance=mol_distance,
                coupling_file=coupling_file,
                x=x,
                y=y,
                Lx=Lx,
                Ly=Ly,
                mol_sites=mol_sites,
            ):
                x2 = (x + Rx) % Lx
                y2 = (y + Ry) % Ly
                s1 = lat(x, y)
                s2 = lat(x2, y2)
                syst[s1, s2] = hop_block

    return syst.finalized()


def report_system_size(fsys: kwant.system.FiniteSystem) -> None:
    H = fsys.hamiltonian_submatrix(sparse=True).tocsr()
    n = H.shape[0]
    nnz = H.nnz
    bytes_total = H.data.nbytes + H.indices.nbytes + H.indptr.nbytes

    print("------ Kwant system report ------")
    print(f"Number of Kwant sites: {len(fsys.sites)}")
    print(f"Orbitals per site: {fsys.sites[0].family.norbs}")
    print(f"Hilbert-space dimension N: {n}")
    print(f"Number of non-zeros (nnz): {nnz}")
    print(f"Sparse H memory: {bytes_total / 1024**2:.2f} MB")
    print(f"Dense H would be: {n * n * 16 / 1024**3:.2f} GB")
    print("---------------------------------")


def estimate_dos(
    fsys: kwant.system.FiniteSystem,
    energy_grid: np.ndarray,
    num_moments: int = 2000,
    rng_seed: int | None = 0,
) -> np.ndarray:
    """Compute DOS via KPM."""
    rng = np.random.default_rng(rng_seed)
    spectrum = kwant.kpm.SpectralDensity(
        fsys,
        num_moments=num_moments,
        rng=rng,
    )
    return spectrum(energy_grid)


if __name__ == "__main__":
    hr = load_hr("wannier90_hr_r0.dat")

    fsys = build_supercell_with_hydrogen(
        hr,
        Lx=10,
        Ly=10,
        coupling_file="hydrogen_interaction_data.h5",
        poly_coeffs_up=[0.0, 0.0, 0.0],
        seed=1234,
    )
    report_system_size(fsys)

    energies = np.linspace(-2.0, 2.0, 400)
    dos = estimate_dos(fsys, energies)
    print("DOS sample:", dos[:5])
