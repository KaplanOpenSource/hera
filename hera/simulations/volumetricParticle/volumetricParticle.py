from decimal import Decimal
import h5py
import warnings
import time
import numpy as np
import os
import pandas as pd
from hera.toolkit import abstractToolkit
from hera import toolkitHome
from hera.datalayer.datahandler import getHandler
import pandas as pd
from calculator import volumetricParticleCalculator
import re

class volumetricParticleAnalysis:
    """
    Analysis layer for the VolumetricParticle toolkit.
    Provides methods for particle data extraction, grid generation (uniform and adaptive),
    and computation of volumetric metrics like concentration, second moment, and occupied volume.
    """

    def __init__(self, toolkit):
        """
        Initialize the analysis layer.

        Args:
            toolkit (VolumetricParticle): The parent toolkit instance.
        """
        self.toolkit = toolkit
        self.calculator = volumetricParticleCalculator()

    def _getParticles(self, t, i):
        """
        Retrieve particle data for a given time step.
        Dispatches to the appropriate source-specific method (e.g., getParticles_openFOAM).

        Args:
            t (float): The time step.

        Returns:
            tuple: A tuple containing (active_ids, positions, mass, age, sigma, TL, epsilon).
        """
        return self.__getattribute__(f"_getParticles_{self.toolkit.source}")(t=t, timestep = i)

    def _getParticles_openFOAM(self, t, **kwargs):
        """
        Retrieve particle data from an OpenFOAM simulation via the OpenFOAM toolkit.

        Args:
            t (float): The time step.

        Returns:
            tuple: (ids, positions, mass, age, sigma, TL, epsilon) as numpy arrays.
        """
        return np.load(os.path.join(self.toolkit.path,f"{float(t):g}",f"origId.npy")).astype(int).ravel(), \
               np.load(os.path.join(self.toolkit.path,f"{float(t):g}",f"globalPositions.npy")), \
               np.load(os.path.join(self.toolkit.path,f"{float(t):g}",f"mass.npy")).ravel(), \
               np.load(os.path.join(self.toolkit.path,f"{float(t):g}",f"age.npy")).ravel(), \
               np.load(os.path.join(self.toolkit.path,f"{float(t):g}",f"sigma.npy")).ravel(), \
               np.load(os.path.join(self.toolkit.path,f"{float(t):g}",f"TL.npy")).ravel(), \
               np.load(os.path.join(self.toolkit.path,f"{float(t):g}",f"epsilon.npy")).ravel()

    def _getParticles_FORTRAN(self, timestep, **kwargs):
        """
        Retrieve particle data from a FORTRAN simulation.

        Args:
            t (float): The time step.

        Returns:
            tuple: (ids, positions, mass, age, sigma, TL, epsilon) as numpy arrays.
        """
        data = pd.read_csv(os.path.join(self.toolkit.path, f"tozaot/machsan/OUTP3d03_3_{timestep:05d}"),
                           names = ["id","mass","age","x","y","z"], delim_whitespace=True)
        data_turb = pd.read_csv(os.path.join(self.toolkit.path, f"tozaot/machsan/OUTTurbProp03_3_{timestep:05d}"),
                           names = ["sigma","TL","epsilon"], delim_whitespace=True)
        return data['id'].values.astype(int), data[['x', 'y', 'z']].values, \
               data['mass'].values, data['age'].values, \
               data_turb['sigma'].values, data_turb['TL'].values, \
               data_turb['epsilon'].values

    def get_uniform_grid(self, domain_edges, dxdydz, overwrite=False, **kwargs):
        """
        Resolves a uniform grid resource. Searches the database for a matching 
        fingerprint; if not found (or overwrite=True), computes and caches a new one.
        
        Args:
            domain_edges (list): [xmin, xmax, ymin, ymax, zmin, zmax].
            dxdydz (float|tuple): Voxel size (isotropic scalar or dx, dy, dz).
            overwrite (bool): If True, forces recalculation and DB update.
            **kwargs: Additional metadata for the fingerprint (e.g., source_id).
            
        Returns:
            CacheDocument: The document handle for the uniform grid.
        """
        # 1. Define the fingerprint for the database search
        fingerprint = {
            "domain_edges": domain_edges,
            "dxdydz": dxdydz,
            **kwargs
        }

        # 2. Try to retrieve from Database
        if not overwrite:
            docs = self.toolkit.getCacheDocuments(type=self.toolkit.UNIFORM_GRID, **fingerprint)
            if docs:
                return docs[0]

        # 3. Computation Logic (Merged from create_uniform_grid)
        xmin, xmax, ymin, ymax, zmin, zmax = domain_edges
        if isinstance(dxdydz, (int, float)):
            dx = dy = dz = float(dxdydz)
        else:
            dx, dy, dz = dxdydz

        # Use a fixed resolution by aligning to multiples of dx/dy/dz and extending beyond maxes
        ox = np.floor(xmin / dx) * dx
        nx = int(np.ceil((xmax - ox) / dx))
        x_edges = ox + np.arange(nx + 1) * dx

        oy = np.floor(ymin / dy) * dy
        ny = int(np.ceil((ymax - oy) / dy))
        y_edges = oy + np.arange(ny + 1) * dy

        oz = np.floor(zmin / dz) * dz
        nz = int(np.ceil((zmax - oz) / dz))
        z_edges = oz + np.arange(nz + 1) * dz

        grid_data = [x_edges, y_edges, z_edges]

        # 4. Resolve unique path and save via Handler
        default_path = os.path.join(self.toolkit.path, self.toolkit.uniformGridFileName)
        filepath = default_path if overwrite else self._get_unique_path(default_path)

        handler = getHandler(self.toolkit.UNIFORM_GRID)
        handler.saveData(resource=grid_data, fileName=filepath)

        # 5. Register and return the document
        doc = self.toolkit.addCacheDocument(
            type=self.toolkit.UNIFORM_GRID,
            resource=filepath,
            dataFormat=self.toolkit.UNIFORM_GRID,
            desc=fingerprint
        )
        return doc

    def get_fixed_res_adaptive_grid(self, df=None, dxdydz=0.01, n_blocks_x=10, domain_edges=None,
            buffer=0.10, overwrite=False, **kwargs):
        """
        Adaptive grid with fixed resolution, where only the Y and Z bounds adapt per X-block.
        """
        if isinstance(dxdydz, (int, float)):
            dx = dy = dz = float(dxdydz)
        else:
            dx, dy, dz = dxdydz

        fingerprint = {
            "dxdydz": dxdydz,
            "n_blocks_x": n_blocks_x,
            "buffer": buffer,
            **kwargs
        }
        default_path = os.path.join(self.toolkit.path, self.toolkit.adaptiveGridFileName)
        filepath = default_path if overwrite else self.toolkit._get_unique_path(default_path)

        docs = self.toolkit.getCacheDocuments(type=self.toolkit.ADAPTIVE_GRID, **fingerprint)
        if docs and not overwrite:
            return docs

        if df is None:
            raise ValueError("Dataframe 'df' must be provided to compute the adaptive grid.")

        if domain_edges is None:
            mins = df[['x', 'y', 'z']].min().compute()
            maxs = df[['x', 'y', 'z']].max().compute()
            global_x_min, global_x_max = mins['x'], maxs['x']
            global_y_min, global_y_max = mins['y'], maxs['y']
            global_z_min, global_z_max = mins['z'], maxs['z']
        else:
            global_x_min, global_x_max, global_y_min, global_y_max, global_z_min, global_z_max = domain_edges

        global_origin = np.floor(np.array([global_x_min, global_y_min, global_z_min]) / dx) * dx

        # Define a global X grid with resolution dx, extending beyond xmax if needed
        nx = int(np.ceil((global_x_max - global_origin[0]) / dx))
        global_x_edges = global_origin[0] + np.arange(nx + 1) * dx

        # Indices to split the nx cells into n_blocks_x blocks
        if n_blocks_x > nx:
            warnings.warn(f"n_blocks_x ({n_blocks_x}) is larger than the total number of cells nx ({nx}). Reducing n_blocks_x to {nx}.")
            n_blocks_x = nx
        split_indices = np.linspace(0, nx, n_blocks_x + 1).astype(int)
        block_boundaries = global_x_edges[split_indices]

        # Assign particles to blocks based on the block boundaries
        def assign_block(pdf):
            # Clip to ensure all particles fall within [0, n_blocks_x - 1]
            pdf['block_id'] = np.clip(np.digitize(pdf['x'], block_boundaries) - 1, 0, n_blocks_x - 1)
            return pdf

        df_with_id = df.map_partitions(assign_block)
        block_stats = df_with_id.groupby('block_id')[['y', 'z']].agg(['min', 'max']).compute()

        master_results = []
        for idx in range(n_blocks_x):
            # X edges for this block are the subset of the global X grid
            start_idx = split_indices[idx]
            end_idx = split_indices[idx+1]
            block_x_edges = global_x_edges[start_idx : end_idx + 1]

            # Ensure the block has at least one cell (should be true if nx >= n_blocks_x)
            if len(block_x_edges) < 2:
                # In the unlikely case of a block with no cells, we create a dummy range
                # to avoid histogramdd failing, though this block will have no particles.
                block_x_edges = np.array([block_boundaries[idx], block_boundaries[idx+1]])

            if idx not in block_stats.index:
                b_y_min, b_y_max = global_y_min, global_y_max
                b_z_min, b_z_max = global_z_min, global_z_max
            else:
                row = block_stats.loc[idx]
                b_y_min, b_y_max = row[('y', 'min')], row[('y', 'max')]
                b_z_min, b_z_max = row[('z', 'min')], row[('z', 'max')]

            # Y alignment
            y_span = max(b_y_max - b_y_min, dy)
            y_center = 0.5 * (b_y_min + b_y_max)
            buffered_y_span = y_span * (1.0 + buffer)
            y_min_buf, y_max_buf = y_center - 0.5 * buffered_y_span, y_center + 0.5 * buffered_y_span
            y_start_aligned = global_origin[1] + np.floor((y_min_buf - global_origin[1]) / dy) * dy
            num_y_cells = max(1, int(np.ceil((y_max_buf - y_start_aligned) / dy)))
            y_edges = y_start_aligned + np.arange(num_y_cells + 1) * dy

            # Z alignment
            z_span = max(b_z_max - b_z_min, dz)
            z_center = 0.5 * (b_z_min + b_z_max)
            buffered_z_span = z_span * (1.0 + buffer)
            z_min_buf, z_max_buf = z_center - 0.5 * buffered_z_span, z_center + 0.5 * buffered_z_span
            z_start_aligned = global_origin[2] + np.floor((z_min_buf - global_origin[2]) / dz) * dz
            num_z_cells = max(1, int(np.ceil((z_max_buf - z_start_aligned) / dz)))
            z_edges = z_start_aligned + np.arange(num_z_cells + 1) * dz

            if domain_edges is not None:
                y_edges = np.unique(np.clip(y_edges, global_y_min, global_y_max))
                z_edges = np.unique(np.clip(z_edges, global_z_min, global_z_max))
                if len(y_edges) < 2: y_edges = np.array([global_y_min, global_y_min + dy])
                if len(z_edges) < 2: z_edges = np.array([global_z_min, global_z_min + dz])

            master_results.append({
                'x_edges': block_x_edges,
                'y_edges': y_edges,
                'z_edges': z_edges,
                'voxel_size': dx # voxel_size usually implies isotropic, here just a marker
            })

        handler = getHandler(self.toolkit.ADAPTIVE_GRID)
        handler.saveData(resource=master_results, fileName=filepath)
        if not docs:
            docs = [self.toolkit.addCacheDocument(
                type=self.toolkit.ADAPTIVE_GRID,
                resource=filepath,
                dataFormat=self.toolkit.ADAPTIVE_GRID,
                desc=fingerprint
            )]
        return docs

    def _scaleExceedingVolumes(self, positions, vol_active, occupied_volumes, grid, block_indices=None, valid_mask=None):
        """
        Scales down particle volumes if the occupied volume fraction in their
        respective grid cell exceeds 1.0. Also clamps the occupied volume grid to [0, 1].
        """
        is_uniform = (isinstance(grid, (list, tuple)) and len(grid) == 3 and not isinstance(grid[0], dict))
        if is_uniform:
            x_edges, y_edges, z_edges = grid
            ix, iy, iz = np.digitize(positions[:, 0], x_edges)-1, \
                np.digitize(positions[:, 1], y_edges)-1, np.digitize(positions[:, 2], z_edges)-1
            valid = (ix >= 0) & (ix < occupied_volumes.shape[0]) & (iy >= 0) & \
                    (iy < occupied_volumes.shape[1]) & (iz >= 0) & (iz < occupied_volumes.shape[2])
            particle_occ = np.ones(len(positions))
            particle_occ[valid] = occupied_volumes[ix[valid], iy[valid], iz[valid]]
            exceed_mask = valid & (particle_occ > 1.0)
            if np.any(exceed_mask):
                vol_active[exceed_mask] *= (1.0 / particle_occ[exceed_mask])
            return np.minimum(occupied_volumes, 1.0)
        else:
            # Use pre-calculated block assignment if available
            if block_indices is None or valid_mask is None:
                x_bounds = np.array([b['x_edges'][0] for b in grid] + [grid[-1]['x_edges'][-1]])
                block_indices = np.digitize(positions[:, 0], x_bounds) - 1
                block_indices[positions[:, 0] == x_bounds[-1]] = len(grid) - 1
                valid_mask = (block_indices >= 0) & (block_indices < len(grid))

            clamped_occ_list = []
            for slice_idx, block_grid in enumerate(occupied_volumes):
                block_info = grid[slice_idx]
                x_edges, y_edges, z_edges = block_info['x_edges'], block_info['y_edges'], block_info['z_edges']

                in_block = (block_indices == slice_idx) & valid_mask
                if not np.any(in_block):
                    clamped_occ_list.append(np.minimum(block_grid, 1.0))
                    continue

                block_pos = positions[in_block]
                # Use searchsorted instead of digitize for faster local mapping
                ix = np.searchsorted(x_edges, block_pos[:, 0], side='right') - 1
                iy = np.searchsorted(y_edges, block_pos[:, 1], side='right') - 1
                iz = np.searchsorted(z_edges, block_pos[:, 2], side='right') - 1

                valid_mask_local = (ix >= 0) & (ix < block_grid.shape[0]) & \
                                  (iy >= 0) & (iy < block_grid.shape[1]) & \
                                  (iz >= 0) & (iz < block_grid.shape[2])

                if np.any(valid_mask_local):
                    p_occ = block_grid[ix[valid_mask_local], iy[valid_mask_local], iz[valid_mask_local]]
                    exceed_mask_local = p_occ > 1.0
                    if np.any(exceed_mask_local):
                        # Map local indices back to global vol_active
                        all_indices = np.where(in_block)[0]
                        valid_indices = all_indices[valid_mask_local]
                        exceed_indices = valid_indices[exceed_mask_local]
                        vol_active[exceed_indices] *= (1.0 / p_occ[exceed_mask_local])

                clamped_occ_list.append(np.minimum(block_grid, 1.0))
            return clamped_occ_list
    
    def preProcessCase(self, nTimesteps):

        self.__getattribute__(f"_preProcessCase_{self.toolkit.source}")(nTimesteps)

    def _preProcessCase_openFOAM(self, nTimesteps):

        def read_foam_field(filename):
            with open(filename, 'r') as f:
                lines = f.readlines()

            for i, line in enumerate(lines):
                match = re.fullmatch(r'\s*(\d+)(?:\s*\{(.*)\})?\s*', line)

                if match:
                    N = int(match.group(1))
                    value = match.group(2)

                    if value is not None:
                        # Uniform field
                        value = value.strip().strip('()')
                        values = np.fromstring(value, sep=' ')
                        return np.tile(values, (N, 1)) if len(values) > 1 else np.full(N, values[0])

                    # Explicit field
                    data = lines[i+1:]
                    if data[0].strip() == '(':
                        data = data[1:]

                    return np.fromstring(' '.join(data[:N]).replace('(', '').replace(')', ''), sep=' ').reshape(N, -1)

        for i in range(1, nTimesteps+1):
            t = Decimal(format(i * self.toolkit.dt, '.12g')).normalize()
            for field in ["origId","globalPositions","mass","age","TL","sigma","epsilon"]:
                np.save(os.path.join(self.toolkit.path,f"{float(t):g}",f"{field}.npy"), \
                        read_foam_field(os.path.join(self.toolkit.path, f"{float(t):g}", "lagrangian", "kinematicCloud", field)))
            print(f"Processed time step {i}")

    def _initialize_accumulator(self, reference):
        arrays = [
            np.zeros(reference.shape, dtype=np.float64)
            for _ in range(2)
        ]
        return arrays

    def _accumulate(self, acc, C, C2):
        sum_C, sum_C2 = acc
        np.add(sum_C, C, out=sum_C)
        np.add(sum_C2, C2, out=sum_C2)


    def _save_interval(self, f, t_str, n_steps, is_uniform, occupied_blocks, block_accumulators):

        if n_steps == 0:
            return

        if is_uniform:
            if uniform_accumulator is None:
                return
            grp_t = f.require_group(t_str)
            for name, arr in zip(self.toolkit.dataset_names, uniform_accumulator):
                if name in grp_t:
                    del grp_t[name]

                # Divide in place to avoid allocating another array.
                np.divide(arr, n_steps, out=arr)
                grp_t.create_dataset(name, data=arr, compression="lzf")

                # Reset for the next interval.
                arr.fill(0)
        else:
            if not occupied_blocks:
                return

            grp_t = f.require_group(t_str)

            for idx in sorted(occupied_blocks):
                acc = block_accumulators[idx]
                slice_name = f"slice_{idx:03d}"

                if slice_name in grp_t:
                    del grp_t[slice_name]
                grp_slice = grp_t.create_group(slice_name)

                for name, arr in zip(self.toolkit.dataset_names, acc):
                    np.divide(arr, n_steps, out=arr)
                    grp_slice.create_dataset(name, data=arr, compression="lzf")
                    arr.fill(0)

            occupied_blocks.clear()

    def calculate(self, grid_doc, nParticles, nTimesteps, initialVolume, sigma0,
                Cr=0.3, mu=0.47, Cphi=1, savedt=1, fileName=None, tm_mode='mean',
                printEvery=100, overwrite=False, debugMode=False, **kwargs):
        """
        Calculate volumetric particle metrics and save temporal averages.
        Saves averages every `savedt` timesteps:
            concentration
            second_moment

        For adaptive grids, unoccupied blocks contribute zero to the
        temporal averages.

        A final incomplete averaging interval is also saved.
        """

        existing_C = self.toolkit.getCacheDocuments(type="concentrations", grid_params=grid_doc.desc, nTimesteps=nTimesteps)
        if not existing_C:
            grid = grid_doc.getData()
            is_uniform = (isinstance(grid, (list, tuple)) and len(grid) == 3 and not isinstance(grid[0], dict))
            results_file = os.path.join(self.toolkit.path, self.toolkit.CONCENTRATION_CACHE)
            results_file = self.toolkit._get_unique_path(results_file)
            with h5py.File(results_file, "w") as f:
                for i in range(1, nTimesteps + 1):

                    t = Decimal(format(i * self.toolkit.dt, ".12g")).normalize()
                    t_str = f"{float(t):g}"

                    # --- Particle loading ---

                    (
                        active_ids,
                        positions,
                        mass,
                        age,
                        sigma_sqrd,
                        particlesTL,
                        particlesEpsilon
                    ) = self._getParticles(t, i)

                    # Empty timesteps still count towards savedt.
                    if len(active_ids) > 0:

                        block_indices = None
                        valid_mask = None

                        if not is_uniform:

                            x_bounds = np.array([b["x_edges"][0] for b in grid] + [grid[-1]["x_edges"][-1]])

                            block_indices = np.digitize(positions[:, 0], x_bounds) - 1
                            block_indices[positions[:, 0] == x_bounds[-1]] = len(grid) - 1
                            valid_mask = ((block_indices >= 0) & (block_indices < len(grid)))

                        results, phi_mean = (
                            self.calculator.calculateConcentration(
                                positions,
                                mass,
                                grid,
                                block_indices=block_indices,
                                valid_mask=valid_mask
                            )
                        )

                        # Save directly to h5 file to avoid keeping all timesteps in memory
                        grp = f.create_group(t_str)
                        grp.create_dataset("phi_mean", data=phi_mean, compression="lzf")
                        if is_uniform:
                            grp.create_dataset("grid", data=results, compression="lzf")
                        else:
                            grid_grp = grp.create_group("grid")
                            for idx, res in enumerate(results):
                                grid_grp.create_dataset(f"slice_{idx:03d}", data=res, compression="lzf")

                    if i % printEvery == 0:
                        print(f"Completed concentration calculation timestep {i} (t={t})")

            self.toolkit.addCacheDocument(
                    type="concentrations",
                    resource=results_file,
                    dataFormat="h5",
                    desc=dict(grid_params=grid_doc.desc, nTimesteps=nTimesteps)
                )

        concentration_cache_file = existing_C[0].resource if existing_C else results_file

        return self.__getattribute__(f"_calculate_{self.toolkit.model}")(
                grid_doc=grid_doc, nParticles=nParticles, nTimesteps=nTimesteps, initialVolume=initialVolume, sigma0=sigma0,
                Cr=Cr, mu=mu, Cphi=Cphi, savedt=savedt, fileName=fileName, tm_mode=tm_mode,
                printEvery=printEvery, overwrite=overwrite, debugMode=debugMode, concentration_cache=concentration_cache_file, **kwargs)

    def _calculate_VPA(self, grid_doc, nParticles, nTimesteps, initialVolume, sigma0,
                Cr=0.3, mu=0.47, Cphi=1, savedt=1,
                scale_exceeding_volumes=True, fileName=None, tm_mode='mean',
                printEvery=100, overwrite=False, debugMode=False, concentration_cache=None, **kwargs):

        # ---------------------------------------------------------
        # 1. Initialization
        # ---------------------------------------------------------

        if not isinstance(savedt, (int, np.integer)) or savedt < 1:
            raise ValueError("savedt must be a positive integer")

        grid = grid_doc.getData()
        is_uniform = (isinstance(grid, (list, tuple)) and len(grid) == 3 and not isinstance(grid[0], dict))
        if fileName is None:
            fileName = (self.toolkit.fileName["uniform"] if is_uniform else self.toolkit.fileName["adaptive"])
        results_file = os.path.join(self.toolkit.path, fileName)
        results_file = self.toolkit._get_unique_path(results_file)

        sim_fingerprint = {
            "grid_params": grid_doc.desc,
            "physics": {
                "dt": self.toolkit.dt,
                "initialVolume": initialVolume,
                "sigma0": sigma0,
                "Cr": Cr,
                "mu": mu,
                "Cphi": Cphi,
            },
            "config": {
                "nParticles": nParticles,
                "nTimesteps": nTimesteps,
                "savedt": savedt,
                "scale_exceeding_volumes": scale_exceeding_volumes,
                "tm_mode": tm_mode,
                "model": "VPA"
            },
            **kwargs
        }

        existing_sims = self.toolkit.getCacheDocuments(type="VolumetricSimulation", **sim_fingerprint)
        if existing_sims:
            if overwrite:
                results_file = existing_sims[0]["desc"]["config"]["results_file"]
            else:
                print(f"Found existing document: {existing_sims[0].resource}")
                return existing_sims

        sim_fingerprint["config"]["results_file"] = results_file

        sigma0_sqrd = sigma0**2
        volume = np.ones(nParticles) * initialVolume
        d_r_sqrd = np.ones(nParticles) * sigma0_sqrd
        sigma_r_sqrd = np.ones(nParticles) * sigma0_sqrd
        t0 = np.zeros(nParticles)

        # ---------------------------------------------------------
        # 2. Accumulator initialization
        # ---------------------------------------------------------

        uniform_accumulator = None

        # Adaptive grids: allocate only when a block first becomes occupied.
        block_accumulators = {} if not is_uniform else None

        # Blocks that received particles during the current averaging interval.
        # This avoids scanning all accumulator arrays when saving.
        occupied_blocks = set()

        interval_steps = 0

        # ---------------------------------------------------------
        # 3. Main computation loop
        # ---------------------------------------------------------

        with h5py.File(results_file, "w") as f:

            for i in range(1, nTimesteps + 1):

                t = Decimal(format(i * self.toolkit.dt, ".12g")).normalize()
                t_str = f"{float(t):g}"

                # --- Particle loading ---

                (
                    active_ids,
                    positions,
                    mass,
                    age,
                    sigma_sqrd,
                    particlesTL,
                    particlesEpsilon
                ) = self._getParticles(t, i)

                # Empty timesteps still count towards savedt.
                if len(active_ids) > 0:

                    # --- State extraction ---

                    vol_active = volume[active_ids]
                    d_r_active = d_r_sqrd[active_ids]
                    sig_r_active = sigma_r_sqrd[active_ids]
                    t0_active = t0[active_ids]

                    # --- Adaptive block assignment ---

                    block_indices = None
                    valid_mask = None

                    if not is_uniform:

                        x_bounds = np.array([b["x_edges"][0] for b in grid] + [grid[-1]["x_edges"][-1]])

                        block_indices = np.digitize(positions[:, 0], x_bounds) - 1
                        block_indices[positions[:, 0] == x_bounds[-1]] = len(grid) - 1
                        valid_mask = ((block_indices >= 0) & (block_indices < len(grid)))

                    # --- STEP 1: Concentration ---

                    results, phi_mean = getHandler("concentrations").getData(resource=concentration_cache, timestep=t_str, is_uniform=is_uniform)

                    # --- STEP 2: Physics and volume evolution ---

                    t0_active, d_r_active, sig_r_active, tm = self.calculator.updatePhysicsState(
                        t0_active=t0_active,
                        d_r_active=d_r_active,
                        sig_r_active=sig_r_active,
                        sigma0_sqrd=sigma0_sqrd,
                        Cr=Cr,
                        mu=mu,
                        Cphi=Cphi,
                        dt=self.toolkit.dt,
                        particlesEpsilon=particlesEpsilon,
                        age=age,
                        sigma_sqrd=sigma_sqrd,
                        particlesTL=particlesTL,
                        warn=True
                    )

                    tm = self.calculator.adjustTmByCell(
                        positions, tm, grid,
                        block_indices=block_indices,
                        valid_mask=valid_mask,
                        mode=tm_mode
                    )

                    phi = mass / vol_active
                    delta_phi = (-self.toolkit.dt * (phi - phi_mean) / tm)

                    vol_active = (vol_active * phi / (phi + delta_phi))

                    # --- STEP 3: Occupied volume and scaling ---

                    if scale_exceeding_volumes:
                        occupied_volumes = (
                            self.calculator.calculateOccupiedVolume(
                                positions,
                                vol_active,
                                grid,
                                block_indices=block_indices,
                                valid_mask=valid_mask
                            )
                        )
                        occupied_volumes = self._scaleExceedingVolumes(
                            positions,
                            vol_active,
                            occupied_volumes,
                            grid,
                            block_indices=block_indices,
                            valid_mask=valid_mask
                        )

                    # --- STEP 4: Second moment ---

                    second_moments = (
                        self.calculator.calculateSecondMoment(
                            positions,
                            mass,
                            vol_active,
                            grid,
                            block_indices=block_indices,
                            valid_mask=valid_mask
                        )
                    )

                    # --- Update particle states ---

                    volume[active_ids] = vol_active
                    d_r_sqrd[active_ids] = d_r_active
                    sigma_r_sqrd[active_ids] = sig_r_active
                    t0[active_ids] = t0_active

                    # -------------------------------------------------
                    # STEP 5: Accumulate instantaneous fields and metrics
                    # -------------------------------------------------

                    if is_uniform:
                        if uniform_accumulator is None:
                            uniform_accumulator = (self._initialize_accumulator(results))

                        self._accumulate(uniform_accumulator, results, second_moments)
                    else:
                        for idx, (C, C2) in enumerate(zip(results, second_moments)):

                            # An unoccupied block contributes zero.
                            # Its accumulator is deliberately left untouched.
                            if not np.any(C):
                                continue
                            if idx not in block_accumulators:
                                block_accumulators[idx] = (self._initialize_accumulator(C))
                            self._accumulate(block_accumulators[idx], C, C2)
                            occupied_blocks.add(idx)

                    # --- Debug output ---

                    if debugMode:
                        timestep_dir = os.path.join(self.toolkit.path, t_str)
                        os.makedirs(timestep_dir, exist_ok=True)
                        np.save(os.path.join(timestep_dir, "sig_r.npy"),sig_r_active**0.5)
                        np.save(os.path.join(timestep_dir, "tm.npy"), tm)

                # -------------------------------------------------
                # STEP 6: Interval management and saving
                # -------------------------------------------------

                interval_steps += 1
                if interval_steps == savedt:
                    self._save_interval(f, t_str, interval_steps, is_uniform, occupied_blocks, block_accumulators)
                    interval_steps = 0
                if i % printEvery == 0:
                    print(f"Completed timestep {i} (t={t})")

            # ---------------------------------------------------------
            # 4. Save final incomplete interval
            # ---------------------------------------------------------

            if interval_steps > 0:
                self._save_interval(f, t_str, interval_steps, is_uniform, occupied_blocks, block_accumulators)

        # ---------------------------------------------------------
        # 5. Finalize
        # ---------------------------------------------------------

        if not existing_sims:
            existing_sims = [
                self.toolkit.addCacheDocument(
                    type="VolumetricSimulation",
                    resource=results_file,
                    dataFormat="h5",
                    desc=sim_fingerprint
                )
            ]

        return existing_sims

    def _extract_profile_data(self, master_grid, results, profile_dim, fixed_dim1, val1, fixed_dim2, val2):
        """
        Internal helper to extract 1D profile data.
        profile_dim: 'x', 'y', or 'z'
        fixed_dim1/2: the other two dimensions
        val1/2: value or [low, high] range
        """
        def get_range(val):
            is_range = isinstance(val, (list, tuple, np.ndarray)) and len(val) == 2
            low, high = (val[0], val[1]) if is_range else (val, val)
            if low > high: low, high = high, low
            return low, high, is_range

        low1, high1, is_range1 = get_range(val1)
        low2, high2, is_range2 = get_range(val2)

        is_list_format = isinstance(results, list)
        has_extra = False
        if is_list_format:
            has_extra = 'occupied_volume' in results[0] if results else False
        elif isinstance(results, dict):
            has_extra = 'occupied_volume' in results and results['occupied_volume'] is not None

        raw_profiles = []
        edges_list = []

        for idx, slice_info in enumerate(master_grid):
            if not isinstance(slice_info, dict): continue
            edges_map = {
                'x': slice_info.get('x_edges'),
                'y': slice_info.get('y_edges'),
                'z': slice_info.get('z_edges')
            }
            if any(v is None for v in edges_map.values()): continue

            if is_list_format:
                res_dict = results[idx]
                concentration = res_dict.get('concentration')
                second_moment = res_dict.get('second_moment')
            else:
                concentration = results.get('concentration')[idx] if results.get('concentration') else None
                second_moment = results.get('second_moment')[idx] if results.get('second_moment') else None

            if concentration is None: continue
            conc_grid = concentration
            sec_grid = second_moment

            # Find matching indices for fixed dimensions
            def find_matching(dim_name, l, h, is_rng):
                edges = edges_map[dim_name]
                matches = []
                for i in range(len(edges) - 1):
                    if not is_rng:
                        if edges[i] <= l < edges[i+1] or (i == len(edges)-2 and l == edges[-1]):
                            matches.append(i)
                    elif not (edges[i+1] <= l or edges[i] >= h):
                        matches.append(i)
                return matches

            match1 = find_matching(fixed_dim1, low1, high1, is_range1)
            match2 = find_matching(fixed_dim2, low2, high2, is_range2)
            if not match1 or not match2: continue

            # Resolve profile dimension
            prof_edges = edges_map[profile_dim]
            n_prof = len(prof_edges) - 1
            s_conc = np.full(n_prof, np.nan)
            s_sec = np.full(n_prof, np.nan)

            for ip in range(n_prof):
                vals_c, vals_s= [], []
                for i1 in match1:
                    for i2 in match2:
                        # Map coordinates based on profile_dim
                        coords = {'x': 0, 'y': 1, 'z': 2}
                        idx_map = [0, 0, 0]
                        idx_map[coords[profile_dim]] = ip
                        idx_map[coords[fixed_dim1]] = i1
                        idx_map[coords[fixed_dim2]] = i2

                        c_idx = tuple(idx_map)
                        if 0 <= c_idx[0] < conc_grid.shape[0] and 0 <= c_idx[1] < conc_grid.shape[1] and 0 <= c_idx[2] < conc_grid.shape[2]:
                            vals_c.append(conc_grid[c_idx])
                            if sec_grid is not None and 0 <= c_idx[0] < sec_grid.shape[0] and 0 <= c_idx[1] < sec_grid.shape[1] and 0 <= c_idx[2] < sec_grid.shape[2]:
                                vals_s.append(sec_grid[c_idx])

                if vals_c:
                    s_conc[ip] = np.nanmean(vals_c)
                    if sec_grid is not None: s_sec[ip] = np.nanmean(vals_s) if vals_s else np.nan

            raw_profiles.append({'edges': prof_edges, 'conc': s_conc, 'sec': s_sec})
            edges_list.append(prof_edges)

        if not raw_profiles:
            return (np.array([]), np.array([]))

        # Create a global reference grid spanning the entire domain at the coarsest resolution
        all_edges_flat = np.concatenate(edges_list)
        global_min, global_max = np.min(all_edges_flat), np.max(all_edges_flat)
        coarsest_res = max([np.mean(np.diff(e)) for e in edges_list])

        ref_edges = np.arange(global_min, global_max + coarsest_res, coarsest_res)
        ref_edges[-1] = global_max
        ref_centers = 0.5 * (ref_edges[:-1] + ref_edges[1:])

        collected_c, collected_s = [], []
        for prof in raw_profiles:
            e_p = prof['edges']
            if np.array_equal(e_p, ref_edges):
                collected_c.append(prof['conc']); collected_s.append(prof['sec'])
            else:
                c_res, s_res = np.full_like(ref_centers, np.nan), np.full_like(ref_centers, np.nan)
                centers_p, dp = 0.5 * (e_p[:-1] + e_p[1:]), np.diff(e_p)
                for i_ref in range(len(ref_edges) - 1):
                    mi, ma = ref_edges[i_ref], ref_edges[i_ref + 1]
                    mask = (centers_p >= mi) & (centers_p < ma)
                    if np.any(mask):
                        w = dp[mask]
                        c_res[i_ref] = np.sum(prof['conc'][mask] * w) / np.sum(w)
                        s_res[i_ref] = np.sum(prof['sec'][mask] * w) / np.sum(w)

                collected_c.append(c_res); collected_s.append(s_res)

        avg_c, avg_s = np.nanmean(collected_c, axis=0), np.nanmean(collected_s, axis=0)

        return ref_centers, avg_c, avg_s

    def get_profiles_x(self, master_grid, results, y0=0.0, z0=0.0):
        """
        Extracts and returns the concentration, second moment, and optional occupied volume
        profiles along the x-axis at a specific (y0, z0) coordinate or averaged over ranges (y0, z0).
        """
        return self._extract_profile_data(master_grid, results, 'x', 'y', y0, 'z', z0)

    def get_profiles_y(self, master_grid, results, x0=0.0, z0=0.0):
        """
        Extracts and returns the y-coordinate profile along the y-axis.
        """
        return self._extract_profile_data(master_grid, results, 'y', 'x', x0, 'z', z0)

    def get_profiles_z(self, master_grid, results, x0=0.0, y0=0.0):
        """
        Extracts vertical profiles along z at (x0, y0) or ranges.
        """
        return self._extract_profile_data(master_grid, results, 'z', 'x', x0, 'y', y0)

    def _extract_slice_data(self, master_grid, results, fixed_dim, fixed_val):
        """
        Internal helper to extract 2D slice data.
        fixed_dim: 'x', 'y', or 'z'
        fixed_val: value or [low, high] range
        """

        # 1. Setup range for the fixed dimension
        is_range = isinstance(fixed_val, (list, tuple, np.ndarray)) and len(fixed_val) == 2
        low, high = (fixed_val[0], fixed_val[1]) if is_range else (fixed_val, fixed_val)
        if low > high: low, high = high, low

        is_list_format = isinstance(results, list)

        # We will store results as a list of DataFrames (one per block) and then concatenate
        block_dfs = []

        for idx, slice_info in enumerate(master_grid):
            if not isinstance(slice_info, dict): continue
            x_edges = slice_info.get('x_edges')
            y_edges = slice_info.get('y_edges')
            z_edges = slice_info.get('z_edges')
            if x_edges is None or y_edges is None or z_edges is None: continue

            # Extract data for this block
            if is_list_format:
                res_dict = results[idx]
                concentration = res_dict.get('concentration')
                second_moment = res_dict.get('second_moment')
            else:
                concentration = results.get('concentration')[idx] if results.get('concentration') else None
                second_moment = results.get('second_moment')[idx] if results.get('second_moment') else None

            if concentration is None: continue
            conc_grid = concentration
            sec_grid = second_moment

            # Find matching indices for the fixed dimension
            edges = {'x': x_edges, 'y': y_edges, 'z': z_edges}[fixed_dim]
            matching_indices = []
            for i in range(len(edges) - 1):
                if not is_range:
                    if edges[i] <= low < edges[i+1] or (i == len(edges)-2 and low == edges[-1]):
                        matching_indices.append(i)
                elif not (edges[i+1] <= low or edges[i] >= high):
                    matching_indices.append(i)

            if not matching_indices: continue

            # Resolve the other two dimensions
            if fixed_dim == 'x':
                dim1_edges, dim2_edges = y_edges, z_edges
                dim1_name, dim2_name = 'y', 'z'
            elif fixed_dim == 'y':
                dim1_edges, dim2_edges = x_edges, z_edges
                dim1_name, dim2_name = 'x', 'z'
            else: # fixed_dim == 'z'
                dim1_edges, dim2_edges = x_edges, y_edges
                dim1_name, dim2_name = 'x', 'y'

            # Create the 2D arrays for this block
            ny, nz = len(dim1_edges)-1, len(dim2_edges)-1 # Generic names for the other 2 dims
            s_conc, s_sec = np.zeros((ny, nz)), np.zeros((ny, nz)) if sec_grid is not None else None

            for i1 in range(ny):
                for i2 in range(nz):
                    vals_c, vals_s = [], []
                    # Loop over the fixed dimension indices
                    for ifix in matching_indices:
                        # Map 3D coords based on fixed_dim
                        if fixed_dim == 'x':
                            coords = (ifix, i1, i2)
                        elif fixed_dim == 'y':
                            coords = (i1, ifix, i2)
                        else:
                            coords = (i1, i2, ifix)
                        if 0 <= coords[0] < conc_grid.shape[0] and 0 <= coords[1] < conc_grid.shape[1] and 0 <= coords[2] < conc_grid.shape[2]:
                            vals_c.append(conc_grid[coords])
                            if sec_grid is not None and 0 <= coords[0] < sec_grid.shape[0] and 0 <= coords[1] < sec_grid.shape[1] and 0 <= coords[2] < sec_grid.shape[2]:
                                vals_s.append(sec_grid[coords])

                    if vals_c:
                        s_conc[i1, i2] = np.nanmean(vals_c)
                        if s_sec is not None: s_sec[i1, i2] = np.nanmean(vals_s) if vals_s else np.nan
                    else:
                        s_conc[i1, i2] = np.nan
                        if s_sec is not None: s_sec[i1, i2] = np.nan

            # Convert this block's 2D slice to a DataFrame
            d1_centers = 0.5 * (dim1_edges[:-1] + dim1_edges[1:])
            d2_centers = 0.5 * (dim2_edges[:-1] + dim2_edges[1:])

            block_data = []
            for i1 in range(ny):
                for i2 in range(nz):
                    row = {dim1_name: d1_centers[i1], dim2_name: d2_centers[i2]}
                    row['concentration'] = s_conc[i1, i2]
                    row['second_moment'] = s_sec[i1, i2] if s_sec is not None else np.nan
                    block_data.append(row)

            block_dfs.append(pd.DataFrame(block_data))

        return pd.concat(block_dfs, ignore_index=True)

    def get_slice_xy(self, master_grid, results, z0=0.0):
        """
        Returns a 2D slice in the XY plane at z0 (or averaged over range z0).
        """
        return self._extract_slice_data(master_grid, results, 'z', z0)

    def get_slice_xz(self, master_grid, results, y0=0.0):
        """
        Returns a 2D slice in the XZ plane at y0 (or averaged over range y0).
        """
        return self._extract_slice_data(master_grid, results, 'y', y0)

    def get_slice_yz(self, master_grid, results, x0=0.0):
        """
        Returns a 2D slice in the YZ plane at x0 (or averaged over range x0).
        """
        return self._extract_slice_data(master_grid, results, 'x', x0)

    def get_max_cell_grids_over_time(self, sim_doc=None, grid_doc=None, h5_file=None, timesteps_list=None,
                                    printEvery=100, overwrite=False, **kwargs):
        """
        Resolves the maximum cell grids over time for an adaptive-grid simulation.

        Compatible with simulations saving temporal averages of:
            concentration
            second_moment

        Returns the maximum value of each saved field over the requested intervals.
        """

        # ---------------------------------------------------------
        # 1. Resolve resources (sim_doc or grid_doc + h5_file)
        # ---------------------------------------------------------

        if sim_doc is not None:
            results_file = sim_doc.resource
            grid_params = sim_doc.desc.get("grid_params")
        elif h5_file is not None and grid_doc is not None:
            results_file = h5_file
            grid_params = None # We already have the grid_doc
        else:
            raise ValueError("Must provide either 'sim_doc' or both 'grid_doc' and 'h5_file'.")

        # Handle timesteps_list
        if timesteps_list is None:
            raise ValueError("timesteps_list must be provided.")

        # Define fingerprint for the output cache
        # If we have sim_doc, use its desc. Otherwise, we use grid_doc and h5_file info.
        fingerprint = {
            "h5_file": results_file,
            "timesteps": sorted(timesteps_list),
            **kwargs
        }
        if sim_doc:
            fingerprint["sim_params"] = sim_doc.desc
            sim_dir = os.path.dirname(results_file)
        else:
            sim_dir = self.toolkit.path
        default_path = os.path.join(sim_dir, "max_cell_grids.h5")
        filepath = default_path if overwrite else self.toolkit._get_unique_path(default_path)

        docs = self.toolkit.getCacheDocuments(
            type=self.toolkit.MAX_CELL_GRID, **fingerprint
        )

        if docs and not overwrite:
            return docs

        # ---------------------------------------------------------
        # 2. Resolve grid document
        # ---------------------------------------------------------

        if grid_doc is not None:
            target_grid_doc = grid_doc
        elif grid_params is not None:
            grid_docs = self.toolkit.getCacheDocuments(
                type=self.toolkit.ADAPTIVE_GRID, **grid_params
            )
            if not grid_docs:
                # Fallback check for uniform grids just in case
                uniform_docs = self.toolkit.getCacheDocuments(
                    type=self.toolkit.UNIFORM_GRID, **grid_params
                )
                if uniform_docs:
                    raise ValueError("Function is available for adaptive grid only!")
                raise ValueError(
                    f"Could not find a grid document matching the parameters: {grid_params}"
                )
            target_grid_doc = grid_docs[0]
        else:
            raise ValueError("Could not resolve grid document.")

        master_grid = target_grid_doc.getData()

        # ---------------------------------------------------------
        # 3. Initialize maximum arrays
        # ---------------------------------------------------------

        max_results = []

        for slice_info in master_grid:

            nx = len(slice_info["x_edges"]) - 1
            ny = len(slice_info["y_edges"]) - 1
            nz = len(slice_info["z_edges"]) - 1

            max_results.append({
                name: np.zeros((nx, ny, nz), dtype=np.float64)
                for name in self.toolkit.dataset_names
            })

        # ---------------------------------------------------------
        # 4. Read and calculate maxima
        # ---------------------------------------------------------

        count = 0

        with h5py.File(results_file, "r") as f:

            for i in timesteps_list:
                t = Decimal(format(i * self.toolkit.dt, ".12g")).normalize()
                t_str = f"{float(t):g}"

                if t_str in f:
                    grp_t = f[t_str]

                    for idx in range(len(master_grid)):
                        slice_name = f"slice_{idx:03d}"
                        if slice_name not in grp_t:
                            continue

                        grp_slice = grp_t[slice_name]
                        max_slice = max_results[idx]

                        for name in self.toolkit.dataset_names:
                            if name not in grp_slice:
                                continue

                            data = grp_slice[name]

                            values = data[:]
                            np.maximum(max_slice[name], values, out=max_slice[name])

                count += 1

                if count % printEvery== 0:
                    print(f"Finished calculation of maxs for time step number {count}")

        # ---------------------------------------------------------
        # 5. Save and register
        # ---------------------------------------------------------

        handler = getHandler("max_cell_grids")

        handler.saveData(resource=max_results, fileName=filepath)

        if not docs:
            docs = [self.toolkit.addCacheDocument(
                type=self.toolkit.MAX_CELL_GRID,
                resource=filepath,
                dataFormat="max_cell_grids",
                desc=fingerprint
            )]

        return docs

    def get_results_for_timesteps(self, sim_doc=None, h5_file=None, timesteps_list=None):
        """
        Retrieves time-averaged simulation results for a list of saved timesteps.

        Returns a list of dictionaries, with None for timesteps that were not saved.

        Adaptive grids:
            Each field contains a list of arrays, one per block.
            Missing blocks are represented by None.

        Uniform grids:
            Each field contains a single numpy array.
        """

        if sim_doc is not None:
            results_file = sim_doc.resource
        elif h5_file is not None:
            results_file = h5_file
        else:
            raise ValueError("Must provide either 'sim_doc' or 'h5_file'.")

        if timesteps_list is None:
            raise ValueError("timesteps_list must be provided.")

        results_list = []


        with h5py.File(results_file, "r") as f:

            for i in timesteps_list:

                # Match the time string format used in calculate()
                t = Decimal(format(i * self.toolkit.dt, '.12g')).normalize()
                t_str = f"{float(t):g}"
                if t_str not in f:
                    results_list.append(None)
                    continue
                grp = f[t_str]

                # -------------------------------------------------
                # Adaptive grid retrieval
                # -------------------------------------------------

                if any(k.startswith('slice_') for k in grp.keys()):

                    timestep_res = {}
                    n_blocks = grp.attrs.get('n_blocks', 0)
                    slice_keys = sorted(k for k in grp.keys() if k.startswith('slice_'))

                    for field in self.dataset_names:
                        field_data = [None] * n_blocks
                        for sk in slice_keys:
                            idx = int(sk.split('_')[1])
                            if idx < n_blocks:
                                grp_slice = grp[sk]
                                if field in grp_slice:
                                    field_data[idx] = grp_slice[field][:]

                        timestep_res[field] = field_data

                # -------------------------------------------------
                # Uniform grid retrieval
                # -------------------------------------------------

                else:
                    timestep_res = {field: (grp[field][:] if field in grp else None) for field in self.dataset_names}

                results_list.append(timestep_res)

        return results_list

    def _calculate_clipped(self, grid_doc, nParticles, nTimesteps, initialVolume, sigma0,
                        Cr=0.3, mu=0.47, Cphi=1, savedt=1,
                        fileName=None, tm_mode='mean',
                        printEvery=100, overwrite=False, debugMode=False, concentration_cache=None, **kwargs):
        """
        Calculate volumetric particle metrics for the clipped case where volume is constant.
        The particle volume is fixed at `initialVolume`.
        The concentration `phi` is a state variable updated according to:
        phi = phi - dt * (phi - phi_mean) / tm.
        """

        if not isinstance(savedt, (int, np.integer)) or savedt < 1:
            raise ValueError("savedt must be a positive integer")

        grid = grid_doc.getData()
        is_uniform = (isinstance(grid, (list, tuple)) and len(grid) == 3 and not isinstance(grid[0], dict))

        # --- Calculate constant cell volume and f ---
        if is_uniform:
            x_e, y_e, z_e = grid
            cell_volume = (x_e[1] - x_e[0]) * (y_e[1] - y_e[0]) * (z_e[1] - z_e[0])
        else:
            block = grid[0]
            cell_volume = (block['x_edges'][1] - block['x_edges'][0]) * \
                          (block['y_edges'][1] - block['y_edges'][0]) * \
                          (block['z_edges'][1] - block['z_edges'][0])

        particle_prob = initialVolume / cell_volume

        if fileName is None:
            fileName = (self.toolkit.fileName["uniform"] if is_uniform else self.toolkit.fileName["adaptive"])

        results_file = os.path.join(self.toolkit.path, fileName)
        results_file = self.toolkit._get_unique_path(results_file)

        sim_fingerprint = {
            "grid_params": grid_doc.desc,
            "physics": {
                "dt": self.toolkit.dt,
                "initialVolume": initialVolume,
                "sigma0": sigma0,
                "Cr": Cr,
                "mu": mu,
                "Cphi": Cphi,
            },
            "config": {
                "nParticles": nParticles,
                "nTimesteps": nTimesteps,
                "savedt": savedt,
                "tm_mode": tm_mode,
                "model": "clipped"
            },
            **kwargs
        }

        existing_sims = self.toolkit.getCacheDocuments(type="VolumetricSimulation", **sim_fingerprint)

        if existing_sims:
            if overwrite:
                results_file = existing_sims[0]["desc"]["config"]["results_file"]
            else:
                print(f"Found existing document: {existing_sims[0].resource}")
                return existing_sims

        sim_fingerprint["config"]["results_file"] = results_file
        sigma0_sqrd = sigma0**2

        # State variables
        phi = np.zeros(nParticles)
        initialized = np.zeros(nParticles, dtype=bool)
        d_r_sqrd = np.ones(nParticles) * sigma0_sqrd
        sigma_r_sqrd = np.ones(nParticles) * sigma0_sqrd
        t0 = np.zeros(nParticles)

        uniform_accumulator = None
        block_accumulators = {} if not is_uniform else None
        occupied_blocks = set()
        interval_steps = 0

        with h5py.File(results_file, "w") as f:
            for i in range(1, nTimesteps + 1):
                t = Decimal(format(i * self.toolkit.dt, ".12g")).normalize()
                t_str = f"{float(t):g}"

                (
                    active_ids,
                    positions,
                    mass,
                    age,
                    sigma_sqrd,
                    particlesTL,
                    particlesEpsilon
                ) = self._getParticles(t, i)

                if len(active_ids) > 0:
                    # --- Initialize phi for new particles ---
                    newly_active = active_ids[~initialized[active_ids]]
                    if len(newly_active) > 0:
                        mask_new = ~initialized[active_ids]
                        phi[newly_active] = mass[mask_new] / initialVolume
                        initialized[active_ids] = True

                    phi_active = phi[active_ids]
                    d_r_active = d_r_sqrd[active_ids]
                    sig_r_active = sigma_r_sqrd[active_ids]
                    t0_active = t0[active_ids]

                    # --- Adaptive block assignment ---
                    block_indices = None
                    valid_mask = None
                    if not is_uniform:
                        x_bounds = np.array([b["x_edges"][0] for b in grid] + [grid[-1]["x_edges"][-1]])
                        block_indices = np.digitize(positions[:, 0], x_bounds) - 1
                        block_indices[positions[:, 0] == x_bounds[-1]] = len(grid) - 1
                        valid_mask = (block_indices >= 0) & (block_indices < len(grid))

                    # --- STEP 1: Concentration ---
                    results, phi_mean = getHandler("concentrations").getData(resource=concentration_cache, timestep=t_str, is_uniform=is_uniform)

                    # --- STEP 2: Physics and volume evolution ---
                    t0_active, d_r_active, sig_r_active, tm = self.calculator.updatePhysicsState(
                        t0_active=t0_active,
                        d_r_active=d_r_active,
                        sig_r_active=sig_r_active,
                        sigma0_sqrd=sigma0_sqrd,
                        Cr=Cr,
                        mu=mu,
                        Cphi=Cphi,
                        dt=self.toolkit.dt,
                        particlesEpsilon=particlesEpsilon,
                        age=age,
                        sigma_sqrd=sigma_sqrd,
                        particlesTL=particlesTL,
                        warn=False
                    )
                    tm = self.calculator.adjustTmByCell(positions, tm, grid, block_indices=block_indices, valid_mask=valid_mask, mode=tm_mode)

                    # Update phi according to the clipped rule
                    phi_active = phi_active - self.toolkit.dt * (phi_active - phi_mean) / tm
                    phi[active_ids] = phi_active

                    # --- STEP 3: Second moment ---
                    if is_uniform:
                        edges = grid
                        n_particles_grid, _ = np.histogramdd(positions, bins=edges)
                        sum_phi_grid, _ = np.histogramdd(positions, bins=edges, weights=phi_active)
                        sum_phi_sq_grid, _ = np.histogramdd(positions, bins=edges, weights=phi_active**2)

                        ov_grid = (n_particles_grid * initialVolume) / cell_volume
                        # Assuming OV < 1 as per user instruction
                        denom = 1.0 - ov_grid
                        c_low_grid = (results - f * sum_phi_grid) / denom
                        second_moments = (1.0 - ov_grid) * (c_low_grid**2) + f * sum_phi_sq_grid
                    else:
                        second_moments = []
                        for idx, block in enumerate(grid):
                            x_edges, y_edges, z_edges = block['x_edges'], block['y_edges'], block['z_edges']

                            # Mask for particles in this block
                            block_mask = (block_indices == idx) & valid_mask
                            if not np.any(block_mask):
                                second_moments.append(np.zeros((len(x_edges)-1, len(y_edges)-1, len(z_edges)-1)))
                                continue

                            block_pos = positions[block_mask]
                            block_phi = phi_active[block_mask]

                            n_p_block, _ = np.histogramdd(block_pos, bins=[x_edges, y_edges, z_edges])
                            sum_phi_block, _ = np.histogramdd(block_pos, bins=[x_edges, y_edges, z_edges], weights=block_phi)
                            sum_phi_sq_block, _ = np.histogramdd(block_pos, bins=[x_edges, y_edges, z_edges], weights=block_phi**2)

                            ov_block = (n_p_block * initialVolume) / cell_volume
                            # Assuming OV < 1 as per user instruction
                            denom = 1.0 - ov_block

                            # results[idx] is the concentration grid for this block
                            res_grid = np.asarray(results[idx])

                            c_low_block = (res_grid - particle_prob * sum_phi_block) / denom

                            sm_block = (1.0 - ov_block) * (c_low_block**2) + particle_prob * sum_phi_sq_block
                            second_moments.append(sm_block)

                    # Update state
                    d_r_sqrd[active_ids] = d_r_active
                    sigma_r_sqrd[active_ids] = sig_r_active
                    t0[active_ids] = t0_active

                    # --- STEP 4: Accumulate ---
                    if is_uniform:
                        if uniform_accumulator is None:
                            uniform_accumulator = _initialize_accumulator(results)
                        self._accumulate(uniform_accumulator, results, second_moments)
                    else:
                        for idx, (C, C2) in enumerate(zip(results, second_moments)):
                            if not np.any(C): continue
                            if idx not in block_accumulators:
                                block_accumulators[idx] = self._initialize_accumulator(C)
                            self._accumulate(block_accumulators[idx], C, C2)
                            occupied_blocks.add(idx)
                            
                    # --- Debug output ---

                    if debugMode:
                        timestep_dir = os.path.join(self.toolkit.path, t_str)
                        os.makedirs(timestep_dir, exist_ok=True)
                        np.save(os.path.join(timestep_dir, "sig_r.npy"),sig_r_active**0.5)
                        np.save(os.path.join(timestep_dir, "tm.npy"), tm)

                interval_steps += 1
                if interval_steps == savedt:
                    self._save_interval(f, t_str, interval_steps, is_uniform, occupied_blocks, block_accumulators)
                    interval_steps = 0
                if i % printEvery == 0:
                    print(f"Completed timestep {i} (t={t})")

            if interval_steps > 0:
                self._save_interval(f, t_str, interval_steps, is_uniform, occupied_blocks, block_accumulators)

        if not existing_sims:
            existing_sims = [
                self.toolkit.addCacheDocument(
                    type="VolumetricSimulation",
                    resource=results_file,
                    dataFormat="h5",
                    desc=sim_fingerprint
                )
            ]
        return existing_sims

class volumetricParticle(abstractToolkit):
    toolkitName = "VolumetricParticle"

    def __init__(self, projectName, path, dt, source="openFOAM", model="VPA", **kwargs):
        super().__init__(projectName=projectName, toolkitName=self.toolkitName, **kwargs)
        if os.path.isabs(path):
            self.path = path
        else:
            self.path = os.path.join(self.FilesDirectory, path)

        os.makedirs(self.path, exist_ok=True)

        self.dt = dt
        self.source = source
        self.t = dt
        self.of = toolkitHome.getToolkit(toolkitHome.SIMULATIONS_OPENFOAM, projectName=projectName)
        self.fileName = {"uniform":"uniform_grid_results.h5", "adaptive":"adaptive_grid_results.h5"}
        self.adaptiveGridFileName = "adaptive_grid.h5"
        self.uniformGridFileName = "uniform_grid.h5"
        self._analysis = volumetricParticleAnalysis(self)
        self.ADAPTIVE_GRID = "adaptive_grid"
        self.UNIFORM_GRID = "uniform_grid"
        self.MAX_CELL_GRID = "max_cell_grid"
        self.CONCENTRATION_CACHE = "concentration_cache.h5"
        self.dataset_names = ["concentration", "second_moment"]
        self.model = model

    @property
    def analysis(self):
        return self._analysis
    
    def _get_unique_path(self, filepath):
        """
        Checks if a file exists and appends a counter if it does.
        Example: grid.h5 -> grid_1.h5 -> grid_2.h5
        """
        if not os.path.exists(filepath):
            return filepath

        root, ext = os.path.splitext(filepath)
        counter = 1
        while os.path.exists(f"{root}_{counter}{ext}"):
            counter += 1

        return f"{root}_{counter}{ext}"
