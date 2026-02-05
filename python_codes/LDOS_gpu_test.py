import kwant
import kwant.kpm
import cupy as cp
import cupyx.scipy.sparse as cusparse
import numpy as np
from math import pi

class GPUEnhancedKwantKPM:
    """
    Enhances Kwant's KPM with GPU acceleration for the most expensive parts.
    Keeps Kwant's interface while accelerating sparse matrix operations.
    """
    
    def __init__(self, system, operator=None, params=None, 
                 num_moments=100, num_vectors=10, energy_resolution=None,
                 gpu_device=0, **kwargs):
        """
        Initialize with Kwant system, but prepare GPU acceleration.
        
        Parameters match kwant.kpm.SpectralDensity for compatibility.
        """
        # Store for CPU fallback
        self.system = system
        self.operator = operator
        self.params = params or {}
        self.num_moments = num_moments
        self.num_vectors = num_vectors
        
        # Get Hamiltonian as sparse matrix
        self.H_cpu = system.hamiltonian_submatrix(sparse=True, params=params).tocsr()
        self.n = self.H_cpu.shape[0]
        
        # Setup GPU
        self.gpu_device = gpu_device
        self._setup_gpu()
        
        # Create both CPU and GPU KPM objects
        self.cpu_kpm = kwant.kpm.SpectralDensity(
            system, operator=operator, params=params,
            num_moments=num_moments, num_vectors=num_vectors,
            energy_resolution=energy_resolution, **kwargs
        )
        
        print(f"GPU-enhanced KPM initialized on device {gpu_device}")
        print(f"System size: {self.n} orbitals")
        
    def _setup_gpu(self):
        """Setup GPU environment and transfer Hamiltonian."""
        cp.cuda.Device(self.gpu_device).use()
        
        # Transfer Hamiltonian to GPU
        self.H_gpu = cusparse.csr_matrix(self.H_cpu)
        
        # Estimate scaling for Chebyshev expansion
        self.scale = self._estimate_scale()
        self.H_gpu = self.H_gpu / self.scale
        
        # Create CUDA stream for async operations
        self.stream = cp.cuda.Stream()
        
    def _estimate_scale(self):
        """Estimate spectral radius for Chebyshev expansion."""
        # Quick Gershgorin estimate
        diag = self.H_cpu.diagonal()
        offdiag_sum = np.abs(self.H_cpu.sum(axis=1).ravel() - diag)
        radius = np.max(np.abs(diag) + offdiag_sum)
        return radius * 1.1  # 10% safety margin
    
    def _projected_indices(self, operator):
        """
        Extract indices for projection operator.
        For hydrogen PDOS: orbital 0 (up) and 1 (down) at hydrogen sites.
        """
        if operator is None:
            return None, None
            
        # This depends on your operator structure
        # For your hydrogen PDOS: returns indices for up and down orbitals
        # You'll need to adapt this based on your operator function
        
        # Example for your case (adjust as needed):
        norb = self.system.sites[0].family.norbs  # orbitals per site
        
        # Find hydrogen sites (first two orbitals at certain positions)
        # This is simplified - you'll need your actual site identification logic
        hydrogen_sites = []
        for i, site in enumerate(self.system.sites):
            tag = site.tag
            # Check if this is a hydrogen site based on your logic
            # For checkerboard pattern in 10x10 with norbs=506:
            if (tag[0] + tag[1]) % 2 != 0:  # Checkerboard condition
                hydrogen_sites.append(i)
        
        # Hydrogen has orbitals 0 (up) and 1 (down)
        indices_up = [i * norb for i in hydrogen_sites]  # orbital 0
        indices_down = [i * norb + 1 for i in hydrogen_sites]  # orbital 1
        
        return indices_up, indices_down
    
    def gpu_chebyshev_moments(self, indices=None):
        """
        GPU-accelerated Chebyshev moment computation.
        
        Parameters:
        indices: specific orbital indices to project onto (for PDOS)
                 if None, computes full trace (for total DOS)
        
        Returns:
        moments: Chebyshev moments
        """
        n = self.n
        M = self.num_moments
        R = self.num_vectors
        
        if indices is not None:
            # Projected moments
            P = len(indices)
            indices_gpu = cp.array(indices, dtype=cp.int32)
            moments_gpu = cp.zeros((M, P), dtype=cp.complex128)
        else:
            # Full trace
            moments_gpu = cp.zeros(M, dtype=cp.complex128)
        
        print(f"GPU: Computing {M} moments with {R} vectors...")
        
        with self.stream:
            for v in range(R):
                # Random vector on GPU
                phi0 = cp.random.randn(n) + 1j * cp.random.randn(n)
                phi0 = phi0 / cp.linalg.norm(phi0)
                
                if indices is not None:
                    phi0_proj = phi0[indices_gpu]
                
                # Chebyshev recursion
                phi_prev = phi0.copy()
                phi_curr = self.H_gpu.dot(phi0)
                
                # First two moments
                if indices is None:
                    moments_gpu[0] += cp.vdot(phi0, phi_prev)
                    moments_gpu[1] += cp.vdot(phi0, phi_curr)
                else:
                    moments_gpu[0] += cp.conj(phi0_proj) * phi_prev[indices_gpu]
                    moments_gpu[1] += cp.conj(phi0_proj) * phi_curr[indices_gpu]
                
                # Remaining moments
                for m in range(2, M):
                    phi_next = 2.0 * self.H_gpu.dot(phi_curr) - phi_prev
                    
                    if indices is None:
                        moments_gpu[m] += cp.vdot(phi0, phi_next)
                    else:
                        moments_gpu[m] += cp.conj(phi0_proj) * phi_next[indices_gpu]
                    
                    phi_prev, phi_curr = phi_curr, phi_next
        
        self.stream.synchronize()
        
        # Average over random vectors
        moments_gpu = moments_gpu / R
        
        # Apply Jackson kernel
        m = cp.arange(M, dtype=cp.float64)
        g = ((M - m + 1) * cp.cos(pi * m / (M + 1)) + 
             cp.sin(pi * m / (M + 1)) / cp.tan(pi / (M + 1))) / (M + 1)
        
        if indices is None:
            moments_gpu = moments_gpu * g
        else:
            moments_gpu = moments_gpu * g[:, cp.newaxis]
        
        return cp.asnumpy(moments_gpu)
    
    def __call__(self, energy_points=None, broadening=None):
        """
        Compatible interface with Kwant's SpectralDensity.
        Uses GPU acceleration for the heavy computation.
        """
        # Get projected indices based on operator
        indices_up, indices_down = self._projected_indices(self.operator)
        
        if indices_up is None:
            # Total DOS
            moments = self.gpu_chebyshev_moments()
            if energy_points is None:
                energy_points = np.linspace(-self.scale, self.scale, 1000)
            dos = self._moments_to_dos(moments, energy_points, broadening)
            return energy_points, dos
        else:
            # Projected DOS (your hydrogen case)
            # Compute moments for up and down separately or together
            all_indices = indices_up + indices_down
            moments_2d = self.gpu_chebyshev_moments(all_indices)
            
            # Split moments
            P_up = len(indices_up)
            moments_up = moments_2d[:, :P_up]
            moments_down = moments_2d[:, P_up:]
            
            # Get energy points
            if energy_points is None:
                energy_points = np.linspace(-self.scale, self.scale, 1000)
            
            # Compute PDOS for up and down
            pdos_up = self._moments_to_dos(moments_up, energy_points, broadening)
            pdos_down = self._moments_to_dos(moments_down, energy_points, broadening)
            
            return energy_points, pdos_up, pdos_down
    
    def _moments_to_dos(self, moments, energies, broadening=None):
        """Convert moments to DOS."""
        # Scale energies to [-1, 1] interval
        x = energies / self.scale
        x = np.clip(x, -0.9999, 0.9999)
        
        # Use Clenshaw's algorithm for stable evaluation
        if moments.ndim == 1:
            return self._clenshaw_eval(moments, x)
        else:
            # Multiple sets of moments (e.g., for different sites)
            dos = np.zeros((len(energies), moments.shape[1]))
            for i in range(moments.shape[1]):
                dos[:, i] = self._clenshaw_eval(moments[:, i], x)
            return dos
    
    def _clenshaw_eval(self, moments, x):
        """Clenshaw's algorithm for Chebyshev summation."""
        N = len(moments) - 1
        b_kp2 = np.zeros_like(x, dtype=np.complex128)
        b_kp1 = np.zeros_like(x, dtype=np.complex128)
        
        for k in range(N, 1, -1):
            b_k = moments[k] + 2 * x * b_kp1 - b_kp2
            b_kp2, b_kp1 = b_kp1, b_k
        
        result = moments[0] + x * b_kp1 - b_kp2
        
        # Normalization factor for DOS
        return (result.real * 2 / (pi * np.sqrt(1 - x**2 + 1e-12)))
        
def compute_hydrogen_pdos_kpm_gpu_integrated(
    fsys: kwant.system.FiniteSystem,
    mol_sites: dict[tuple[int, int], int],
    *,
    orbital_up: int = 0,
    orbital_down: int = 1,
    energy_grid: np.ndarray | None = None,
    num_moments: int = 1500,
    num_vectors: int = 30,
    gpu_device: int = 0,
    use_kwant_interface: bool = True,
) -> tuple[np.ndarray, list[tuple[int, int]], np.ndarray, np.ndarray]:
    """
    Compute hydrogen PDOS using GPU-accelerated Kwant KPM.
    
    Parameters:
    -----------
    use_kwant_interface: If True, uses Kwant's operator interface.
                         If False, uses direct index projection (faster).
    """
    norb = fsys.sites[0].family.norbs
    coords = [xy for xy, _ in sorted(mol_sites.items(), key=lambda kv: kv[1])]
    num_sites = len(coords)
    
    if use_kwant_interface:
        # Use Kwant's operator interface (more flexible)
        def hydrogen_operator(bra, ket, **params):
            """Operator projecting onto hydrogen orbitals."""
            result_up = []
            result_down = []
            
            for site_idx, site in enumerate(fsys.sites):
                xy = tuple(site.tag)
                if xy in mol_sites:
                    # Hydrogen orbitals 0 (up) and 1 (down)
                    offset = site_idx * norb
                    result_up.append(np.conj(bra[offset + orbital_up]) * ket[offset + orbital_up])
                    result_down.append(np.conj(bra[offset + orbital_down]) * ket[offset + orbital_down])
            
            return np.array(result_up + result_down)
        
        # Create GPU-enhanced KPM with Kwant operator
        spectrum = GPUEnhancedKwantKPM(
            fsys,
            operator=hydrogen_operator,
            num_moments=num_moments,
            num_vectors=num_vectors,
            gpu_device=gpu_device,
        )
        
        # Compute PDOS
        if energy_grid is None:
            energy_grid = np.linspace(-5, 5, 1000)
        
        energies, pdos_up, pdos_down = spectrum(energy_grid)
        
        # Reshape to (energies, sites)
        pdos_up = pdos_up.reshape(len(energies), num_sites)
        pdos_down = pdos_down.reshape(len(energies), num_sites)
        
    else:
        # Direct index projection (faster, less flexible)
        # Build site-to-index mapping
        site_to_id = {site: i for i, site in enumerate(fsys.sites)}
        
        # Collect orbital indices
        indices_up = []
        indices_down = []
        for site in fsys.sites:
            xy = tuple(site.tag)
            if xy in mol_sites:
                sid = site_to_id[site]
                indices_up.append(sid * norb + orbital_up)
                indices_down.append(sid * norb + orbital_down)
        
        # Create KPM without operator (will use direct indices)
        spectrum = GPUEnhancedKwantKPM(
            fsys,
            num_moments=num_moments,
            num_vectors=num_vectors,
            gpu_device=gpu_device,
        )
        
        # Compute moments for all indices at once
        all_indices = indices_up + indices_down
        moments_2d = spectrum.gpu_chebyshev_moments(all_indices)
        
        # Split moments
        moments_up = moments_2d[:, :len(indices_up)]
        moments_down = moments_2d[:, len(indices_up):]
        
        # Compute energy grid
        if energy_grid is None:
            energy_grid = np.linspace(-spectrum.scale, spectrum.scale, 1000)
        
        # Convert moments to PDOS
        pdos_up = spectrum._moments_to_dos(moments_up, energy_grid)
        pdos_down = spectrum._moments_to_dos(moments_down, energy_grid)
        
        energies = energy_grid
    
    return energies, coords, pdos_up, pdos_down

def build_hydrogen_nbp_supercell_optimized(
    hr: HRData,
    Lx: int = 10,
    Ly: int = 10,
    coupling_file: str | Path = "SOC_linregress.h5",
    E_F: float = 5.79371650,
    *,
    use_gpu: bool = False,
    gpu_device: int = 0,
    random_seed: int = None,
) -> tuple[kwant.system.FiniteSystem, dict[tuple[int, int], int]]:
    """
    Optimized builder for hydrogen on NbP supercell.
    
    Parameters:
    -----------
    E_F: Fermi energy for hydrogen onsite energy adjustment
    random_seed: For reproducible hydrogen distances
    """
    import time
    start = time.time()
    
    if random_seed is not None:
        np.random.seed(random_seed)
    
    # Hydrogen parameters
    hydrogen_orbitals = 2
    norb_total = hr.num_wann + hydrogen_orbitals
    lat = kwant.lattice.square(norbs=norb_total)
    syst = kwant.Builder()
    
    # Poly coefficients from your data
    poly_coeffs_up = [-0.5371981017543859, -0.24168142111111088 + E_F]
    poly_coeffs_dn = [-0.46845213684210535, -0.2524031988888888 + E_F]
    
    # Create random distances for hydrogen (-0.6 to 0.6 Å)
    num_hydrogen_sites = (Lx * Ly) // 2  # Checkerboard pattern
    mol_distance = np.random.uniform(-0.6, 0.6, size=hydrogen_orbitals * num_hydrogen_sites)
    
    # Map checkerboard positions to molecule indices
    mol_sites = {}
    mol_idx = 0
    for x in range(Lx):
        for y in range(Ly):
            if (x + y) % 2 != 0:  # Checkerboard condition
                if mol_idx < mol_distance.size:
                    mol_sites[(x, y)] = mol_idx
                    mol_idx += hydrogen_orbitals
    
    print(f"Building {Lx}x{Ly} NbP supercell with hydrogen")
    print(f"  Total sites: {Lx * Ly}")
    print(f"  Hydrogen sites: {len(mol_sites)}")
    print(f"  Total orbitals: {norb_total * Lx * Ly}")
    print(f"  Hydrogen orbitals: {hydrogen_orbitals * len(mol_sites)}")
    
    # Extract onsite Hamiltonian from Wannier90
    onsite_base = np.zeros((norb_total, norb_total), dtype=complex)
    for idx in range(hr.nrpts):
        if np.all(hr.R[:, idx] == 0):
            onsite_base[hydrogen_orbitals:, hydrogen_orbitals:] = hr.H_R[:, :, idx].copy()
            break
    
    # Place onsite terms with hydrogen energy adjustment
    for x in range(Lx):
        for y in range(Ly):
            onsite = onsite_base.copy()
            
            # Adjust hydrogen onsite energies based on distance
            if (x, y) in mol_sites:
                mol_idx = mol_sites[(x, y)]
                onsite[0, 0] = np.polyval(poly_coeffs_up, mol_distance[mol_idx])
                onsite[1, 1] = np.polyval(poly_coeffs_dn, mol_distance[mol_idx + 1])
            
            syst[lat(x, y)] = onsite
    
    # Add hoppings
    print("Adding hoppings...")
    
    # Pre-calculate hydrogen coupling blocks if using GPU
    if use_gpu:
        # Pre-load coupling data to GPU memory
        coupling_data = _preload_coupling_data_gpu(coupling_file, gpu_device)
    else:
        coupling_data = None
    
    for idx in range(hr.nrpts):
        R_vec = tuple(int(v) for v in hr.R[:, idx])
        if R_vec == (0, 0, 0):
            continue
        
        # Get hopping from Wannier90
        hop_base = hr.H_R[:, :, idx] / hr.weight[idx]
        
        # Create full hopping block
        hop_block = np.zeros((norb_total, norb_total), dtype=complex)
        hop_block[hydrogen_orbitals:, hydrogen_orbitals:] = hop_base
        
        # Add hydrogen coupling for relevant sites
        for x in range(Lx):
            for y in range(Ly):
                x2 = (x + R_vec[0]) % Lx
                y2 = (y + R_vec[1]) % Ly
                
                # Check if either site has hydrogen
                src_has_h = (x, y) in mol_sites
                dst_has_h = (x2, y2) in mol_sites
                
                if not src_has_h and not dst_has_h:
                    syst[lat(x, y), lat(x2, y2)] = hop_block
                    continue
                
                # Compute hydrogen coupling
                if use_gpu and coupling_data:
                    hydrogen_coupling = _compute_hydrogen_coupling_gpu(
                        R_vec, x, y, x2, y2, mol_sites, mol_distance,
                        coupling_data, hr.num_wann, hydrogen_orbitals
                    )
                else:
                    hydrogen_coupling = coupling_hydrogen_slab(
                        R_vec, hr.num_wann, mol_distance, coupling_file,
                        x, y, Lx, Ly, mol_sites, use_gpu=False
                    )
                
                syst[lat(x, y), lat(x2, y2)] = hop_block + hydrogen_coupling
    
    print(f"Build time: {time.time() - start:.2f} seconds")
    
    return syst.finalized(), mol_sites

def _preload_coupling_data_gpu(coupling_file, gpu_device):
    """Pre-load hydrogen coupling data to GPU memory."""
    cp.cuda.Device(gpu_device).use()
    
    coupling_data = {}
    with h5py.File(coupling_file, 'r') as f:
        for R_vec_key in f.keys():
            # Parse R vector from key like "[0,1,0]"
            R_vec = tuple(map(int, R_vec_key.strip('[]').split(',')))
            
            block_data = {}
            for label in f[R_vec_key]:
                if label.endswith('_im'):
                    continue
                    
                mapping = coupling_hydrogen_slab_mapping(label)
                if mapping is None:
                    continue
                
                # Load data to GPU
                data_re = cp.asarray(f[R_vec_key][label][:])
                data_im = cp.asarray(f[R_vec_key][f"{label}_im"][:])
                
                block_data[label] = {
                    'mapping': mapping,
                    'data_re': data_re,
                    'data_im': data_im,
                    'indices': cp.asarray(mapping[2])
                }
            
            coupling_data[R_vec] = block_data
    
    return coupling_data
    
if __name__ == "__main__":
    import sys
    import time
    
    # Parse arguments
    use_gpu = "--use-gpu" in sys.argv
    quick_test = "--quick-test" in sys.argv
    args = [arg for arg in sys.argv[1:] if arg not in ["--use-gpu", "--quick-test"]]
    
    if quick_test:
        # Quick test with small system
        hr = load_hr(args[0])
        Lx, Ly = 4, 4  # Smaller for testing
        
        fsys, mol_sites = build_hydrogen_nbp_supercell_optimized(
            hr, Lx, Ly, use_gpu=use_gpu
        )
        
        # Test GPU KPM
        energies, coords, pdos_up, pdos_down = compute_hydrogen_pdos_kpm_gpu_integrated(
            fsys, mol_sites,
            num_moments=500,  # Fewer for testing
            num_vectors=10,
            gpu_device=0,
            use_kwant_interface=False  # Faster for testing
        )
        
        print(f"Quick test complete: {len(energies)} energy points")
        print(f"Hydrogen sites: {len(coords)}")
        
        sys.exit(0)
    
    # Full calculation
    if len(args) < 2:
        print("Usage: python script.py <wannier90_hr.dat> <output.csv> [--use-gpu] [--quick-test]")
        sys.exit(1)
    
    # Load Hamiltonian
    hr = load_hr(args[0])
    
    print("=" * 60)
    print("Hydrogen PDOS on NbP Supercell")
    print("=" * 60)
    print(f"GPU acceleration: {use_gpu}")
    
    start_total = time.time()
    
    # Build system
    fsys, mol_sites = build_hydrogen_nbp_supercell_optimized(
        hr, Lx=10, Ly=10, use_gpu=use_gpu
    )
    
    # Define energy range for PDOS
    E_F = 5.79371650
    energies = np.linspace(-1.06, E_F + 1, 800)
    
    # Compute PDOS
    if use_gpu:
        energies, coords, pdos_up, pdos_down = compute_hydrogen_pdos_kpm_gpu_integrated(
            fsys, mol_sites,
            energy_grid=energies,
            num_moments=1500,
            num_vectors=30,
            gpu_device=0,
            use_kwant_interface=False  # Direct indexing is faster
        )
    else:
        # Fallback to original CPU Kwant KPM
        energies, coords, pdos_up, pdos_down = compute_hydrogen_pdos_kpm(
            fsys, mol_sites,
            energy_grid=energies,
            num_moments=1500,
            num_vectors=30
        )
    
    # Save results
    df = pd.DataFrame({"Energies": energies})
    for idx, (x, y) in enumerate(coords):
        df[f"H_{x}_{y}_up"] = pdos_up[:, idx]
        df[f"H_{x}_{y}_down"] = pdos_down[:, idx]
    
    df.to_csv(args[1], index=False)
    
    print(f"\nTotal execution time: {time.time() - start_total:.2f} seconds")
    print(f"Results saved to {args[1]}")
    print("=" * 60)
