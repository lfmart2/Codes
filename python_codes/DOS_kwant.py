from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import h5py
import kwant
import numpy as np
import sys


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


def _build_hydrogen_coupling_map() -> dict[str, tuple[int, int, np.ndarray]]:
    mapping: dict[str, tuple[int, int, np.ndarray]] = {}
    layer_specs = {
        "first_d_layer": (range(2, 12, 2), range(3, 12, 2)),
        "second_d_layer": (range(12, 22, 2), range(13, 22, 2)),
        "third_d_layer": (range(22, 32, 2), range(23, 32, 2)),
        "first_sp_layer": (range(280, 288, 2), range(281, 288, 2)),
        "second_sp_layer": (range(288, 296, 2), range(289, 296, 2)),
        "third_sp_layer": (range(296, 304, 2), range(297, 304, 2)),
    }
    for prefix, (even_range, odd_range) in layer_specs.items():
        even_idx = np.fromiter(even_range, dtype=int)
        odd_idx = np.fromiter(odd_range, dtype=int)
        for spin_index, spin_label in enumerate(("up", "dn")):
            for target_label, target_idx in (("up", even_idx), ("dn", odd_idx)):
                for suffix in ("", "_im"):
                    label = f"{prefix}_{spin_label}_{target_label}{suffix}"
                    if target_label == "up":
                        mapping[label] = (spin_index, 0, target_idx)
                    elif target_label == "dn":
                        mapping[label] = (spin_index, 1, target_idx)
                    else:
                        raise ValueError(f"Unexpected target label: {target_label}")
    return mapping


HYDROGEN_COUPLING_MAP = _build_hydrogen_coupling_map()


def coupling_hydrogen_slab_mapping(
    label: str,
) -> tuple[int, int, np.ndarray] | None:
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
    hydrogen_orbitals = 2
    final_block = np.zeros(
        (num_slab + hydrogen_orbitals, num_slab + hydrogen_orbitals), dtype=complex
    )
    rows_up = np.array([0], dtype=int)
    rows_dn = rows_up + 1
    src_idx = mol_sites.get((x, y))
    dst_idx = mol_sites.get(((x + R_vec[0]) % Lx, (y + R_vec[1]) % Ly))
    has_src = src_idx is not None
    has_dst = dst_idx is not None
    if not has_src and not has_dst:
        return final_block
    dist_up = (
        np.array([mol_distance[src_idx]], dtype=float) if has_src else np.array([])
    )
    dist_dn = (
        np.array([mol_distance[src_idx + 1]], dtype=float)
        if has_src
        else np.array([])
    )
    dist_up_conj = (
        np.array([mol_distance[dst_idx]], dtype=float) if has_dst else np.array([])
    )
    dist_dn_conj = (
        np.array([mol_distance[dst_idx + 1]], dtype=float)
        if has_dst
        else np.array([])
    )

    list_Ri = _load_hydrogen_block_keys(str(file_path), R_vec)
    with h5py.File(file_path, "r") as file:
        block_interaction = file[f"[{R_vec[0]},{R_vec[1]},{R_vec[2]}]"]
        for label in list_Ri:
            if label.endswith("_im"):
                continue
            mapping = coupling_hydrogen_slab_mapping(label)
            if mapping is None:
                continue
            spin_index, target_label_spin, interaction_indices = mapping
            data_re = block_interaction[label]
            data_im = block_interaction[f"{label}_im"]
            cols = interaction_indices + hydrogen_orbitals
            if has_src:
                if spin_index == 0:
                    distances = dist_up
                    rows = rows_up
                else:
                    distances = dist_dn
                    rows = rows_dn
                values = (
                    distances[:, None] * data_re[:, 1][None, :]
                    + data_re[:, 0][None, :]
                ) + 1j * (
                    distances[:, None] * data_im[:, 1][None, :]
                    + data_im[:, 0][None, :]
                )
                final_block[np.ix_(rows, cols)] += values
            if has_dst:
                if target_label_spin == 0:
                    distances_conj = dist_up_conj
                    rows_conj = rows_up
                else:
                    distances_conj = dist_dn_conj
                    rows_conj = rows_dn
                values_conj = (
                    distances_conj[:, None] * data_re[:, 1][None, :]
                    + data_re[:, 0][None, :]
                ) - 1j * (
                    distances_conj[:, None] * data_im[:, 1][None, :]
                    + data_im[:, 0][None, :]
                )
                final_block[np.ix_(cols, rows_conj)] += values_conj.T
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
    hydrogen_orbitals = 2
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
    poly_coeffs_dn: list[float],
) -> kwant.system.FiniteSystem:
    """Build a 2D supercell including hydrogen coupling data."""
    num_molecules = (Lx * Ly) // 2
    hydrogen_orbitals = 2
    norb = hr.num_wann + hydrogen_orbitals
    lat = kwant.lattice.square(norbs=norb)
    syst = kwant.Builder()

    mol_distance = np.random.uniform(
        -0.6, 0.6, size=hydrogen_orbitals * num_molecules
    )

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
            mol_idx += 2

    for x in range(Lx):
        for y in range(Ly):
            onsite = onsite_block.copy()
            mol_idx = mol_sites.get((x, y))
            if mol_idx is not None:
                onsite[0, 0] = np.polyval(poly_coeffs_up, mol_distance[mol_idx])
                onsite[1, 1] = np.polyval(poly_coeffs_dn, mol_distance[mol_idx + 1])
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


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit("Usage: python DOS_kwant.py <wannier90_hr.dat> <output.csv>")
    hr = load_hr(sys.argv[1])

    coupling_file = "SOC_linregress_by_R.h5"

    poly_coeffs_up = [
        2.0127743055555554,
        -0.2429976190476202,
        -0.09279894841269883,
        5.44317319047619,
    ]
    poly_coeffs_dn = [
        1.5496388888888881,
        -0.19418601190476167,
        0.051376289682539204,
        5.424151619047619,
    ]

    fsys = build_supercell_with_hydrogen(
        hr,
        Lx=10,
        Ly=10,
        coupling_file=coupling_file,
        poly_coeffs_up=poly_coeffs_up,
        poly_coeffs_dn=poly_coeffs_dn,
    )
    report_system_size(fsys)
    rho = kwant.kpm.SpectralDensity(fsys)

    emin, emax = rho.bounds
    Emax_req = 5.53348380 + 1.0
    print(f"Minimum bound energy: {emin}")
    print(f"Minimum bound energy: {emax}")
    energies = np.linspace(emin + 0.01, Emax_req, 40)
    dos = rho(energies)

    filename_np = sys.argv[2]
    np.savetxt(
        filename_np,
        np.column_stack([energies, dos]),
        delimiter=",",
        header="Energies, DOS",
        comments="",
    )
