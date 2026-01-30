from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import h5py
import kwant
import numpy as np
import pandas as pd
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
        # Weight lines may span multiple rows; keep accumulating until nrpts is reached.
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

        # Map each unique lattice vector R -> column index in the arrays.
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

            # Wannier90 indices start at 1; convert to zero-based for the arrays.
            H_R[i - 1, j - 1, idx] = real + 1j * imag

        nr_eff = idx_counter
        return HRData(
            num_wann=num_wann,
            nrpts=nr_eff,
            R=R[:, :nr_eff],
            weight=weight_arr[:nr_eff],
            H_R=H_R[:, :, :nr_eff],
        )
    # update the new Fermi energies accordingly woth the SCF calculations instead of the NSCCF -> Wannier90 calculations


def _build_hydrogen_coupling_map() -> dict[str, tuple[int, int, np.ndarray]]:
    mapping: dict[str, tuple[int, int, np.ndarray]] = {}
    layer_specs = {
        "first_d_layer": (range(0, 10, 2), range(1, 10, 2)),
        "second_d_layer": (range(10, 20, 2), range(11, 20, 2)),
        "third_d_layer": (range(20, 30, 2), range(21, 30, 2)),
        "first_sp_layer": (range(280, 288, 2), range(281, 288, 2)),
        "second_sp_layer": (range(288, 296, 2), range(289, 296, 2)),
        "third_sp_layer": (range(296, 304, 2), range(297, 304, 2)),
    }
    # Build label -> (spin, target_spin, target_indices) mapping once to reuse downstream.
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


# @lru_cache(maxsize=32)
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

    # Distances for polynomial coupling (separate spin channels).
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

            # Source contribution: rows correspond to molecule hydrogen orbitals.
            if has_src:
                if spin_index == 0:
                    distances = dist_up
                    rows = rows_up
                else:
                    distances = dist_dn
                    rows = rows_dn
                values = (
                    data_re[:, 3][None, :] * distances[:, None] ** 3 +
                    data_re[:, 2][None, :] * distances[:, None] ** 2 +
                    data_re[:, 1][None, :] * distances[:, None] + 
                    data_re[:, 0][None, :]
                ) + 1j * (
                    data_im[:, 3][None, :] * distances[:, None] ** 3 +
                    data_im[:, 2][None, :] * distances[:, None] ** 2 +
                    data_im[:, 1][None, :] * distances[:, None] + 
                    data_im[:, 0][None, :]
                )
                final_block[np.ix_(rows, cols)] += values

            # Destination contribution: Hermitian conjugate with target spin index.
            if has_dst:
                if target_label_spin == 0:
                    distances_conj = dist_up_conj
                    rows_conj = rows_up
                else:
                    distances_conj = dist_dn_conj
                    rows_conj = rows_dn
                values_conj = (
                    data_re[:, 3][None, :] * distances_conj[:, None] ** 3 +
                    data_re[:, 2][None, :] * distances_conj[:, None] ** 2 +
                    data_re[:, 1][None, :] * distances_conj[:, None] + 
                    data_re[:, 0][None, :]
                ) + 1j * (
                    data_im[:, 3][None, :] * distances_conj[:, None] ** 3 +
                    data_im[:, 2][None, :] * distances_conj[:, None] ** 2 +
                    data_im[:, 1][None, :] * distances_conj[:, None] + 
                    data_im[:, 0][None, :]
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

        # Base hopping from Wannier90 (normalized by weight).
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

    # Randomize molecule-surface distances in [-0.6, 0.6] Angstrom (uniform).
    mol_distance = np.random.uniform(
        -0.6, 0.6, size=hydrogen_orbitals * num_molecules
    )

    # Extract onsite block (R=0) and copy into the system onsite matrices.
    onsite_block = np.zeros((norb, norb), dtype=complex)
    for idx in range(hr.nrpts):
        if np.all(hr.R[:, idx] == 0):
            onsite_block[hydrogen_orbitals:, hydrogen_orbitals:] = hr.H_R[
                :, :, idx
            ].copy()
            break

    # Map checkerboard positions -> molecule index (two hydrogen orbitals each).
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

    # Place onsite terms and shift hydrogen energies using fitted polynomials.
    for x in range(Lx):
        for y in range(Ly):
            onsite = onsite_block.copy()
            mol_idx = mol_sites.get((x, y))
            if mol_idx is not None:
                onsite[0, 0] = np.polyval(poly_coeffs_up, mol_distance[mol_idx])
                onsite[1, 1] = np.polyval(poly_coeffs_dn, mol_distance[mol_idx + 1])
            syst[lat(x, y)] = onsite

    # Add hoppings within the periodic supercell, including hydrogen coupling.
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

    return syst.finalized(), mol_sites


def compute_hydrogen_pdos_kpm(
    fsys: kwant.system.FiniteSystem,
    mol_sites: dict[tuple[int, int], int],
    *,
    orbital_up: int = 0,
    orbital_down: int = 1,
    energy_grid: np.ndarray | None = None,
    num_moments: int = 1500,
    num_vectors: int = 30,
) -> tuple[np.ndarray, list[tuple[int, int]], np.ndarray, np.ndarray]:
    """Compute hydrogen PDOS for two orbitals using KPM + custom operator.

    Returns
    -------
    energies : (NE,)
    coords : list of (x, y) in mol_sites order
    pdos_up : (NE, Ns)
    pdos_down : (NE, Ns)
    """
    norb = fsys.sites[0].family.norbs
    coords = [xy for xy, _ in sorted(mol_sites.items(), key=lambda kv: kv[1])]
    coords_set = set(coords)

    site_to_id = {site: i for i, site in enumerate(fsys.sites)}
    dof_up: list[int] = []
    dof_down: list[int] = []
    for site in fsys.sites:
        xy = tuple(site.tag)
        if xy in coords_set:
            sid = site_to_id[site]
            dof_up.append(sid * norb + orbital_up)
            dof_down.append(sid * norb + orbital_down)

    dof_up = np.asarray(dof_up, dtype=int)
    dof_down = np.asarray(dof_down, dtype=int)

    if len(dof_up) != len(coords):
        missing = [xy for xy in coords if xy not in {tuple(s.tag) for s in fsys.sites}]
        raise ValueError(
            f"Requested {len(coords)} sites, matched {len(dof_up)}. Missing? {missing[:10]}"
        )

    def op(bra: np.ndarray, ket: np.ndarray, **_params: object) -> np.ndarray:
        v_up = np.conjugate(bra[dof_up]) * ket[dof_up]
        v_down = np.conjugate(bra[dof_down]) * ket[dof_down]
        return np.concatenate([v_up, v_down])

    spectrum = kwant.kpm.SpectralDensity(
        fsys,
        operator=op,
        num_moments=num_moments,
        num_vectors=num_vectors,
        mean=True,
    )
    if energy_grid is None:
        energies, dens = spectrum()
    else:
        energies, dens = spectrum(np.asarray(energy_grid))
    dens = np.asarray(dens)

    num_sites = len(coords)
    dens_2 = dens.reshape(len(energies), 2, num_sites)
    pdos_up = dens_2[:, 0, :]
    pdos_down = dens_2[:, 1, :]

    return energies, coords, pdos_up, pdos_down

if __name__ == "__main__":
    # Expect Wannier90 Hamiltonian path and output filename from CLI.
    if len(sys.argv) < 3:
        raise SystemExit("Usage: python DOS_kwant.py <wannier90_hr.dat> <output.csv>")
    hr = load_hr(sys.argv[1])

    E_F =  5.79371650

    coupling_file = "SOC_linregress.h5"

    # Linear polynomial coefficients for hydrogen onsite shifts (spin-resolved).
    poly_coeffs_up = [
        -0.5371981017543859,
        -0.24168142111111088 + E_F,
    ]
    poly_coeffs_dn = [
        -0.46845213684210535,
        -0.2524031988888888 + E_F,
    ]

    fsys, mol_sites = build_supercell_with_hydrogen(
        hr,
        Lx=10,
        Ly=10,
        coupling_file=coupling_file,
        poly_coeffs_up=poly_coeffs_up,
        poly_coeffs_dn=poly_coeffs_dn,
        )

    filename_np = sys.argv[2]
    emin = -1.06
    emax = E_F + 1
    energies = np.linspace(emin, emax, 800)

    
    energies, coords, pdos_up, pdos_down = compute_hydrogen_pdos_kpm(
        fsys,
        mol_sites,
        energy_grid=energies,
        num_moments=1500,
        num_vectors=30,
    )
    

    df = pd.DataFrame({"Energies": energies})
    for idx, (x, y) in enumerate(coords):
        df[f"H_{x}_{y}_up"] = pdos_up[:, idx]
        df[f"H_{x}_{y}_down"] = pdos_down[:, idx]
    df.to_csv(filename_np, index=False)
