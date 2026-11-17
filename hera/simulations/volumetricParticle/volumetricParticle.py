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
        data = self.toolkit.of.stochasticLagrangian.getCaseResults(
            caseDescriptor=self.toolkit.path,
            withReleaseTimes=True,
            withVelocity=False,
            withTurbProperties=True,
            cloudName="kinematicCloud",
            cache=False,
            timeList=[f"{float(t):g}"],
            useDask=False)

        return data['id'].values.astype(int), data[['x', 'y', 'z']].values, \
               data['mass'].values, data['age'].values, \
               data['sigma'].values, data['TL'].values, \
               data['epsilon'].values

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
        default_path = os.path.join(self.toolkit.path, self.toolkit.adaptiveGridFileName + "_fixed_res.h5")
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

    def _scaleExceedingVolumes(self, positions, vol_active, occupied_volumes, grid):
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
            clamped_occ_list = []
            for slice_idx, block_grid in enumerate(occupied_volumes):
                # Get the corresponding block edges from the grid definition
                block_info = grid[slice_idx]
                x_edges, y_edges, z_edges = block_info['x_edges'], block_info['y_edges'], block_info['z_edges']

                # 1. Identify particles in this block's X-range
                if slice_idx == len(grid) - 1:
                    in_block = (positions[:, 0] >= x_edges[0]) & (positions[:, 0] <= x_edges[-1])
                else:
                    in_block = (positions[:, 0] >= x_edges[0]) & (positions[:, 0] < x_edges[-1])

                if not np.any(in_block):
                    clamped_occ_list.append(np.minimum(block_grid, 1.0))
                    continue

                # 2. Find 3D coordinates for particles in this block
                block_pos = positions[in_block]
                ix = np.digitize(block_pos[:, 0], x_edges) - 1
                iy = np.digitize(block_pos[:, 1], y_edges) - 1
                iz = np.digitize(block_pos[:, 2], z_edges) - 1

                # 3. Filter for particles that actually fall inside the 3D grid bounds
                valid_mask = (ix >= 0) & (ix < block_grid.shape[0]) & \
                            (iy >= 0) & (iy < block_grid.shape[1]) & \
                            (iz >= 0) & (iz < block_grid.shape[2])

                if np.any(valid_mask):
                    # Extract the occupied volume for the valid particles
                    # block_grid is 3D: [X, Y, Z]
                    p_occ = block_grid[ix[valid_mask], iy[valid_mask], iz[valid_mask]]

                    exceed_mask = p_occ > 1.0
                    if np.any(exceed_mask):
                        # Map back to the original vol_active array
                        # we use the in_block mask and then the valid_mask
                        all_indices = np.where(in_block)[0]
                        valid_indices = all_indices[valid_mask]
                        exceed_indices = valid_indices[exceed_mask]

                        # Scale volumes
                        vol_active[exceed_indices] *= (1.0 / p_occ[exceed_mask])

                clamped_occ_list.append(np.minimum(block_grid, 1.0))
            return clamped_occ_list

    def calculate(self, grid_doc, nParticles, nTimesteps, initialVolume, sigma0, Cr = 0.3, mu = 0.47, Cphi = 1,
                  continue_step = None, savedt = 1, scale_exceeding_volumes = True, fileName = None,
                  printEvery=100, overwrite = False, debugMode = False, **kwargs):
        """
        Main computation loop for volumetric particle metrics over a series of timesteps.

        This method implements a physical model for particle volume evolution based on
        Cassiani (2013). It calculates concentration, second moment, and occupied volume
        for each timestep, averaging results between save intervals.

        Args:
            grid_doc (cache document): The grid document (uniform or adaptive).
            nParticles (int): Total number of particles.
            nTimesteps (int): Number of timesteps to simulate.
            initialVolume (float): Initial volume of a particle.
            sigma0 (float): Initial standard deviation of the particle size distribution.
            Cr (float): Coefficient for volume growth (default 0.3).
            mu (float): Growth rate constant (default 0.47).
            Cphi (float): Concentration correction factor (default 1).
            continue_step (int, optional): Step to resume from. If provided, loads state from state.h5.
            savedt (int): Interval (in steps) for saving results to disk (default 1).
            scale_exceeding_volumes (bool): Whether to scale down volumes if occupied volume exceeds 1.0.
            fileName (str, optional): Filename for the output HDF5 file.
            printEvery (int): Frequency of progress printing (default 100).

        Returns:
            None: Results are written directly to HDF5 files in the toolkit path.
        """
        # Construct Simulation Fingerprint
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
                "scale_exceeding_volumes": scale_exceeding_volumes
            },
            **kwargs
        }

        grid = grid_doc.getData()
        sigma0_sqrd = sigma0**2

        # 1. Handle initialization vs resuming from a checkpoint
        if continue_step is not None:
            t_prev = Decimal(format(continue_step * self.toolkit.dt, '.12g')).normalize()
            state_file = os.path.join(self.toolkit.path, str(t_prev), "state.h5")
            with h5py.File(state_file, "r") as f:
                volume, d_r_sqrd, sigma_r_sqrd, t0 = \
                    f['volume'][:], f['d_r_sqrd'][:], f['sigma_r_sqrd'][:], f['t0'][:]
            start_step = continue_step + 1
        else:
            volume, d_r_sqrd, sigma_r_sqrd, t0 = \
                np.ones(nParticles)*initialVolume, np.ones(nParticles)*sigma0_sqrd, \
                np.ones(nParticles)*sigma0_sqrd, np.zeros(nParticles)
            start_step = 1

        # Robust uniform grid check: list/tuple of length 3 where elements are not dicts
        is_uniform = (isinstance(grid, (list, tuple)) and len(grid) == 3 and not isinstance(grid[0], dict))

        if fileName is None:
            fileName = self.toolkit.fileName["uniform"] if is_uniform else self.toolkit.fileName["adaptive"]

        if not overwrite:
            t_first = Decimal(format(start_step * self.toolkit.dt, '.12g')).normalize()
            first_step_dir = os.path.join(self.toolkit.path, f"{float(t_first):g}")
            first_step_file = os.path.join(first_step_dir, fileName)
            unique_path = self.toolkit._get_unique_path(first_step_file)
            fileName = os.path.basename(unique_path)

        # Skip if run is already finished and we aren't overwriting
        existing_sims = self.toolkit.getCacheDocuments(type="VolumetricSimulation", **sim_fingerprint)
        if existing_sims:
            if not overwrite and continue_step is None:
                print(f"Found existing document: {existing_sims[0].resource}")
                return existing_sims

        # Initialize time-averaging accumulators
        sum_concentration, sum_second_moment, sum_occupied_volume, steps_accumulated = None, None, None, 0

        # 2. Loop over timesteps starting from start_step
        for i in range(start_step, nTimesteps+1):
            t = Decimal(format(i * self.toolkit.dt, '.12g')).normalize()
            active_ids, positions, mass, age, sigma_sqrd, particlesTL, particlesEpsilon = self._getParticles(t, i)

            if len(active_ids) == 0: continue

            mass /= nParticles

            # --- VECTORIZED STATE EXTRACTION VIA INDEXING ---
            vol_active, d_r_active, sig_r_active, t0_active = \
                volume[active_ids], d_r_sqrd[active_ids], sigma_r_sqrd[active_ids], t0[active_ids]

            # --- STEP 1: CONCENTRATION MUST BE CALCULATED BEFORE VOLUME UPDATE ---
            results, phi_mean = self.calculator.calculateConcentration(positions, mass, grid)
            zero_mask = (t0_active == 0)

            # --- STEP 2: PHYSICS & VOLUME EVOLUTION UPDATES ---
            t0_active[zero_mask] = (sigma0_sqrd / Cr / particlesEpsilon[zero_mask])**(1./3.)
            d_r_active += 3 * Cr * particlesEpsilon * (t0_active + age)**2 * self.toolkit.dt # Cassiani 2013, Eq. 38
            new_sigma_r_sqrd = d_r_active / (1 + (d_r_active - sigma0_sqrd) / 
                (sigma0_sqrd + 2 * sigma_sqrd * particlesTL * age)) # Cassiani 2013, Eq. 36
            sig_r_active = np.maximum(sig_r_active, new_sigma_r_sqrd)
            L_sqrd = ((1.5 * sigma_sqrd) ** 3) / (particlesEpsilon ** 2)
            sigma_ur_sqrd = np.minimum(sigma_sqrd, sigma_sqrd * (sig_r_active / L_sqrd)**(1./3.)) # Cassiani 2013, Eq. 35
            k = 1.5 * sigma_sqrd
            tm = np.minimum(mu * np.sqrt(sig_r_active / sigma_ur_sqrd), k / particlesEpsilon / Cphi) # Cassiani 2013, Eq. 34
            if tm.min() < self.toolkit.dt:
                print(f"Warning: minimum calculated tm {tm.min()} < dt!")
            tm = np.maximum(tm, self.toolkit.dt + 1e-10) 
            phi = mass / vol_active

            delta_phi = - self.toolkit.dt * (phi - phi_mean) / tm # Cassiani 2013, Eq. 6
            vol_active = vol_active * phi / (phi + delta_phi) # Cassiani 2013, Eq. 7

            # --- STEP 3: OCCUPIED VOLUME & VOLUME SCALING ---
            occupied_volumes = self.calculator.calculateOccupiedVolume(positions, vol_active, grid)
            if scale_exceeding_volumes:
                occupied_volumes = self._scaleExceedingVolumes(positions, vol_active, occupied_volumes, grid)

            # --- STEP 4: SECOND MOMENT (Using scaled vol_active) ---
            second_moments = self.calculator.calculateSecondMoment(positions, mass, vol_active, grid)

            # --- WRITE UPDATED STATES BACK VECTORIALLY ---
            volume[active_ids], d_r_sqrd[active_ids], sigma_r_sqrd[active_ids], t0[active_ids] = \
             vol_active, d_r_active, sig_r_active, t0_active
                
            # --- TIME-AVERAGING ACCUMULATION ---
            if is_uniform:
                if sum_concentration is None: 
                    sum_concentration, sum_second_moment, sum_occupied_volume = \
                        np.zeros_like(results), np.zeros_like(second_moments), np.zeros_like(occupied_volumes)
                sum_concentration += results
                sum_second_moment += second_moments
                sum_occupied_volume += occupied_volumes
            else:
                if sum_concentration is None: sum_concentration, sum_second_moment, sum_occupied_volume = \
                    [np.zeros_like(res) for res in results], [np.zeros_like(sm) for sm in second_moments], \
                    [np.zeros_like(ov) for ov in occupied_volumes]
                for idx, (res_grid, sm_grid, ov_grid) in enumerate(zip(results, second_moments, occupied_volumes)): 
                    sum_concentration[idx] += res_grid
                    sum_second_moment[idx] += sm_grid
                    sum_occupied_volume[idx] += ov_grid

            steps_accumulated += 1

            # --- SAVING & WRITING (Triggered every `savedt` steps) ---
            if (i-1) % savedt == 0 or i == nTimesteps - 1:
                if is_uniform: avg_concentration, avg_second_moment, avg_occupied_volume = \
                    sum_concentration/steps_accumulated, sum_second_moment/steps_accumulated, \
                    sum_occupied_volume/steps_accumulated
                else: avg_concentration, avg_second_moment, avg_occupied_volume = \
                    [sc/steps_accumulated for sc in sum_concentration], \
                    [ss/steps_accumulated for ss in sum_second_moment], \
                    [sov/steps_accumulated for sov in sum_occupied_volume]

                timestep_dir = os.path.join(self.toolkit.path, f"{float(t):g}")
                os.makedirs(timestep_dir, exist_ok=True)
                file_path = os.path.join(timestep_dir, fileName)

                with h5py.File(file_path, "w") as f:
                    if is_uniform:
                        f.create_dataset('concentration', data=avg_concentration, compression="lzf")
                        f.create_dataset('second_moment', data=avg_second_moment, compression="lzf")
                        f.create_dataset('occupied_volume', data=avg_occupied_volume, compression="lzf")
                    else:
                        for idx, (res_grid, sm_grid, ov_grid) in \
                            enumerate(zip(avg_concentration, avg_second_moment, avg_occupied_volume)):
                            grp = f.create_group(f"slice_{idx:03d}")
                            grp.create_dataset('concentration', data=res_grid, compression="lzf")
                            grp.create_dataset('second_moment', data=sm_grid, compression="lzf")
                            grp.create_dataset('occupied_volume', data=ov_grid, compression="lzf")
                with h5py.File(os.path.join(timestep_dir, "state.h5"), "w") as f:
                    f.create_dataset('volume', data=volume, compression="lzf")
                    f.create_dataset('d_r_sqrd', data=d_r_sqrd, compression="lzf")
                    f.create_dataset('sigma_r_sqrd', data=sigma_r_sqrd, compression="lzf")
                    f.create_dataset('t0', data=t0, compression="lzf")

                sum_concentration = sum_second_moment = sum_occupied_volume = None; steps_accumulated = 0
                if debugMode:
                    np.save(os.path.join(timestep_dir, "sig_r.npy"),sig_r_active**0.5)
                    np.save(os.path.join(timestep_dir, "L.npy"),L_sqrd**0.5)
                    np.save(os.path.join(timestep_dir, "tm.npy"),tm)
                    np.save(os.path.join(timestep_dir, "tm_clip.npy"),k / particlesEpsilon / Cphi)
            if i % printEvery == 0: print(f"Completed timestep {i} (t={t})")

        sim_fingerprint["config"]["filename"] = fileName
        sim_fingerprint["config"]["sim_path"] = self.toolkit.path
        if not existing_sims:
            existing_sims = [self.toolkit.addCacheDocument(
                type="VolumetricSimulation",
                resource=self.toolkit.path,
                dataFormat="folder",
                desc=sim_fingerprint)]

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
        has_occupied_volume = False
        if is_list_format:
            has_occupied_volume = 'occupied_volume' in results[0] if results else False
        elif isinstance(results, dict):
            has_occupied_volume = 'occupied_volume' in results and results['occupied_volume'] is not None

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
                occupied_volume = res_dict.get('occupied_volume')
            else:
                concentration = results.get('concentration')[idx] if results.get('concentration') else None
                second_moment = results.get('second_moment')[idx] if results.get('second_moment') else None
                occupied_volume = results.get('occupied_volume')[idx] if results.get('occupied_volume') else None

            if concentration is None: continue
            conc_grid = concentration
            sec_grid = second_moment
            occ_grid = occupied_volume

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
            s_occ = np.full(n_prof, np.nan)

            for ip in range(n_prof):
                vals_c, vals_s, vals_o = [], [], []
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
                            if occ_grid is not None and 0 <= c_idx[0] < occ_grid.shape[0] and 0 <= c_idx[1] < occ_grid.shape[1] and 0 <= c_idx[2] < occ_grid.shape[2]:
                                vals_o.append(occ_grid[c_idx])

                if vals_c:
                    s_conc[ip] = np.nanmean(vals_c)
                    if sec_grid is not None: s_sec[ip] = np.nanmean(vals_s) if vals_s else np.nan
                    if occ_grid is not None: s_occ[ip] = np.nanmean(vals_o) if vals_o else np.nan

            raw_profiles.append({'edges': prof_edges, 'conc': s_conc, 'sec': s_sec, 'occ': s_occ})
            edges_list.append(prof_edges)

        if not raw_profiles:
            return (np.array([]), np.array([]), np.array([]), np.array([])) if has_occupied_volume else (np.array([]), np.array([]), np.array([]))

        # Create a global reference grid spanning the entire domain at the coarsest resolution
        all_edges_flat = np.concatenate(edges_list)
        global_min, global_max = np.min(all_edges_flat), np.max(all_edges_flat)
        coarsest_res = max([np.mean(np.diff(e)) for e in edges_list])

        ref_edges = np.arange(global_min, global_max + coarsest_res, coarsest_res)
        ref_edges[-1] = global_max
        ref_centers = 0.5 * (ref_edges[:-1] + ref_edges[1:])

        collected_c, collected_s, collected_o = [], [], []
        for prof in raw_profiles:
            e_p = prof['edges']
            if np.array_equal(e_p, ref_edges):
                collected_c.append(prof['conc']); collected_s.append(prof['sec']); collected_o.append(prof['occ'])
            else:
                c_res, s_res, o_res = np.full_like(ref_centers, np.nan), np.full_like(ref_centers, np.nan), np.full_like(ref_centers, np.nan)
                centers_p, dp = 0.5 * (e_p[:-1] + e_p[1:]), np.diff(e_p)
                for i_ref in range(len(ref_edges) - 1):
                    mi, ma = ref_edges[i_ref], ref_edges[i_ref + 1]
                    mask = (centers_p >= mi) & (centers_p < ma)
                    if np.any(mask):
                        w = dp[mask]
                        c_res[i_ref] = np.sum(prof['conc'][mask] * w) / np.sum(w)
                        s_res[i_ref] = np.sum(prof['sec'][mask] * w) / np.sum(w)
                        o_res[i_ref] = np.sum(prof['occ'][mask] * w) / np.sum(w)
                collected_c.append(c_res); collected_s.append(s_res); collected_o.append(o_res)

        avg_c, avg_s, avg_o = np.nanmean(collected_c, axis=0), np.nanmean(collected_s, axis=0), np.nanmean(collected_o, axis=0)
        if has_occupied_volume: return ref_centers, avg_c, avg_s, avg_o
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
        has_occupied_volume = False
        if is_list_format:
            has_occupied_volume = 'occupied_volume' in results[0] if results else False
        elif isinstance(results, dict):
            has_occupied_volume = 'occupied_volume' in results and results['occupied_volume'] is not None

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
                occupied_volume = res_dict.get('occupied_volume')
            else:
                concentration = results.get('concentration')[idx] if results.get('concentration') else None
                second_moment = results.get('second_moment')[idx] if results.get('second_moment') else None
                occupied_volume = results.get('occupied_volume')[idx] if results.get('occupied_volume') else None

            if concentration is None: continue
            conc_grid = concentration
            sec_grid = second_moment
            occ_grid = occupied_volume

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
            s_conc, s_sec, s_occ = np.zeros((ny, nz)), np.zeros((ny, nz)) if sec_grid is not None else None, np.zeros((ny, nz)) if occ_grid is not None else None

            for i1 in range(ny):
                for i2 in range(nz):
                    vals_c, vals_s, vals_o = [], [], []
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
                            if occ_grid is not None and 0 <= coords[0] < occ_grid.shape[0] and 0 <= coords[1] < occ_grid.shape[1] and 0 <= coords[2] < occ_grid.shape[2]:
                                vals_o.append(occ_grid[coords])

                    if vals_c:
                        s_conc[i1, i2] = np.nanmean(vals_c)
                        if s_sec is not None: s_sec[i1, i2] = np.nanmean(vals_s) if vals_s else np.nan
                        if s_occ is not None: s_occ[i1, i2] = np.nanmean(vals_o) if vals_o else np.nan
                    else:
                        s_conc[i1, i2] = np.nan
                        if s_sec is not None: s_sec[i1, i2] = np.nan
                        if s_occ is not None: s_occ[i1, i2] = np.nan

            # Convert this block's 2D slice to a DataFrame
            d1_centers = 0.5 * (dim1_edges[:-1] + dim1_edges[1:])
            d2_centers = 0.5 * (dim2_edges[:-1] + dim2_edges[1:])

            block_data = []
            for i1 in range(ny):
                for i2 in range(nz):
                    row = {dim1_name: d1_centers[i1], dim2_name: d2_centers[i2]}
                    row['concentration'] = s_conc[i1, i2]
                    row['second_moment'] = s_sec[i1, i2] if s_sec is not None else np.nan
                    row['occupied_volume'] = s_occ[i1, i2] if s_occ is not None else np.nan
                    block_data.append(row)

            block_dfs.append(pd.DataFrame(block_data))

        if not block_dfs:
            return pd.DataFrame(columns=[dim1_name, dim2_name, 'concentration', 'second_moment', 'occupied_volume'])

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

    def get_max_cell_grids_over_time(self, sim_doc, timesteps_list, printRate=None, overwrite=False, **kwargs):
        """
        Resolves the max cell grids resource for a specific simulation.
        """
        # 1. Define the fingerprint
        # Use sim_doc.resource as the unique identifier for this simulation run
        fingerprint = {
            "sim_params": sim_doc.desc,
            "timesteps": sorted(timesteps_list),
            **kwargs
        }

        # 2. Try to retrieve from Database
        sim_path = sim_doc.resource
        default_path = os.path.join(sim_path, "max_cell_grids.h5")
        filepath = default_path if overwrite else self.toolkit._get_unique_path(default_path)

        docs = self.toolkit.getCacheDocuments(type=self.toolkit.MAX_CELL_GRID, **fingerprint)
        if docs:
            if not overwrite:
                return docs

        # 3. Resolve the Grid document using the grid_params stored in the sim_doc
        grid_params = sim_doc.desc.get("grid_params")
        if not grid_params:
            raise ValueError("The simulation document is missing 'grid_params' in its description.")

        # Find the grid document that matches the parameters used for this simulation
        grid_docs = self.toolkit.getCacheDocuments(type=self.toolkit.ADAPTIVE_GRID, **grid_params)
        if not grid_docs:
            grid_docs = self.toolkit.getCacheDocuments(type=self.toolkit.UNIFORM_GRID, **grid_params)
            if not grid_docs:
                raise ValueError(f"Could not find a grid document matching the parameters: {grid_params}")
            else:
                raise ValueError("Function is available for adaptive grid only!")

        grid_doc = grid_docs[0]
        master_grid = grid_doc.getData()

        # 4. Computation Logic
        
        fileName = sim_doc.desc["config"]["filename"]

        n_slices = len(master_grid)
        max_results = []
        for slice_info in master_grid:
            # Now we include nx for the 3D block
            nx = len(slice_info['x_edges']) - 1
            ny = len(slice_info['y_edges']) - 1
            nz = len(slice_info['z_edges']) - 1
            max_results.append({
                'concentration': np.zeros((nx, ny, nz)),
                'second_moment': np.zeros((nx, ny, nz)),
                'occupied_volume': np.zeros((nx, ny, nz))
            })

        count = 0
        for i in timesteps_list:
            t = Decimal(format(i * self.toolkit.dt, '.12g')).normalize()
            result_file = os.path.join(sim_path, f"{float(t):g}", fileName)

            if not os.path.exists(result_file):
                continue

            with h5py.File(result_file, "r") as f:
                for idx in range(n_slices):
                    slice_name = f"slice_{idx:03d}"
                    if slice_name in f:
                        grp = f[slice_name]
                        if 'concentration' in grp:
                            c_grid = grp['concentration'][:]
                            if c_grid.size > 0:
                                max_results[idx]['concentration'] = np.maximum(max_results[idx]['concentration'], c_grid)
                        if 'second_moment' in grp:
                            c2_grid = grp['second_moment'][:]
                            if c2_grid.size > 0:
                                max_results[idx]['second_moment'] = np.maximum(max_results[idx]['second_moment'], c2_grid)
                        if 'occupied_volume' in grp:
                            ov_grid = grp['occupied_volume'][:]
                            if ov_grid.size > 0:
                                max_results[idx]['occupied_volume'] = np.maximum(max_results[idx]['occupied_volume'], ov_grid)

            if printRate is not None and count % printRate == 0:
                print(f"Finished calculation for time step number {count}")
            count += 1
        if not docs:
            docs = [self.toolkit.addCacheDocument(
                type=self.toolkit.MAX_CELL_GRID,
                resource=filepath,
                dataFormat="max_cell_grids",
                desc=fingerprint
            )]

        handler = getHandler("max_cell_grids")
        handler.saveData(resource=max_results, fileName=filepath)

        return docs

    def getResultsForTimesteps(self, sim_doc, timesteps_list):
        """
        Retrieves simulation results (concentration, second moment, occupied volume) 
        for a list of timesteps from the simulation's HDF5 files.

        Args:
            sim_doc (CacheDocument): The simulation cache document.
            timesteps_list (list): List of timestep indices (integers).

        Returns:
            list: A list of dictionaries containing the results for each timestep. 
                Returns None for timesteps where the file is missing.
        """
        sim_path = sim_doc.resource
        fileName = sim_doc.desc["config"]["filename"]
        results_list = []

        for i in timesteps_list:
            # Match the exact time string format used during the 'calculate' loop
            t = Decimal(format(i * self.toolkit.dt, '.12g')).normalize()
            result_file = os.path.join(sim_path, str(t), fileName)

            if not os.path.exists(result_file):
                results_list.append(None)
                continue

            with h5py.File(result_file, "r") as f:
                # Determine if the result is adaptive (contains slice groups) or uniform
                if "slice_000" in f:
                    # --- Adaptive Grid Retrieval ---
                    timestep_res = {}
                    # Identify all available slices in sorted order
                    slice_keys = sorted([k for k in f.keys() if k.startswith('slice_')])

                    for field in ['concentration', 'second_moment', 'occupied_volume']:
                        field_data = []
                        for sk in slice_keys:
                            grp = f[sk]
                            field_data.append(grp[field][:] if field in grp else None)
                        timestep_res[field] = field_data
                else:
                    # --- Uniform Grid Retrieval ---
                    timestep_res = {}
                    for field in ['concentration', 'second_moment', 'occupied_volume']:
                        timestep_res[field] = f[field][:] if field in f else None

                results_list.append(timestep_res)

        return results_list

class volumetricParticle(abstractToolkit):
    toolkitName = "VolumetricParticle"

    def __init__(self, projectName, path, dt, source="openFOAM", **kwargs):
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
        self.adaptiveGridFileName = "adaptive_grid"
        self.uniformGridFileName = "uniform_grid.h5"
        self._analysis = volumetricParticleAnalysis(self)
        self.ADAPTIVE_GRID = "adaptive_grid"
        self.UNIFORM_GRID = "uniform_grid"
        self.MAX_CELL_GRID = "max_cell_grid"

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
