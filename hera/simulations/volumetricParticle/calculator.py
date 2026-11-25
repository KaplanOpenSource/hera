import numpy as np

class volumetricParticleCalculator:


    def calculateConcentration(self, original_positions, mass, grid, block_indices=None, valid_mask=None):
        """
        Calculate the mass concentration across the grid.
        Dispatches to either uniform or non-uniform calculation based on the grid type.

        Args:
            original_positions (np.ndarray): Array of particle positions [x, y, z].
            mass (np.ndarray): Array of particle masses.
            grid (list|tuple|list[dict]): The grid definition.
            block_indices (np.ndarray, optional): Pre-calculated block assignment.
            valid_mask (np.ndarray, optional): Pre-calculated validity mask.
        """
        is_uniform = (isinstance(grid, (list, tuple)) and len(grid) == 3 and not isinstance(grid[0], dict))
        return self.calculateConcentration_uniform(original_positions, mass, grid) \
            if is_uniform else self.calculateConcentration_nonuniform(original_positions, mass, grid, block_indices, valid_mask)

    def calculateConcentration_uniform(self, original_positions, mass, grid):
        """
        Calculate concentration for a uniform grid using histogramdd.
        """
        original_positions = np.asarray(original_positions)
        mass = np.asarray(mass).ravel()

        # Unpack grid edges
        x_edges, y_edges, z_edges = grid
        nx = len(x_edges) - 1
        ny = len(y_edges) - 1
        nz = len(z_edges) - 1

        # 1. 3D Histogram with weights (mass) across all dimensions simultaneously
        mass_grid, _ = np.histogramdd(
            original_positions,
            bins=[x_edges, y_edges, z_edges],
            weights=mass
        )

        # 2. Compute 3D cell volumes for rectilinear grids (supports non-cubic or stretched grids)
        dx = np.diff(x_edges)
        dy = np.diff(y_edges)
        dz = np.diff(z_edges)
        cell_volumes = dx[:, None, None] * dy[None, :, None] * dz[None, None, :]

        # 3. Convert mass to concentration (Mass / Cell Volume)
        concentration_grid = mass_grid / cell_volumes

        # 4. Look up concentrations for each particle in original order
        px = original_positions[:, 0]
        py = original_positions[:, 1]
        pz = original_positions[:, 2]

        ix = np.searchsorted(x_edges, px, side='right') - 1
        iy = np.searchsorted(y_edges, py, side='right') - 1
        iz = np.searchsorted(z_edges, pz, side='right') - 1

        # Check which particles fall inside the grid bounds
        valid_mask = (
            (ix >= 0) & (ix < nx) &
            (iy >= 0) & (iy < ny) &
            (iz >= 0) & (iz < nz)
        )

        original_order_concentrations = np.zeros_like(mass)

        if np.any(valid_mask):
            v_idx = np.where(valid_mask)[0]
            original_order_concentrations[v_idx] = concentration_grid[
                ix[v_idx], iy[v_idx], iz[v_idx]
            ]

        return concentration_grid, original_order_concentrations

    def calculateConcentration_nonuniform(self, original_positions, mass, master_grid, block_indices=None, valid_mask=None):
        """
        Calculate concentration for a non-uniform (adaptive) grid using 3D blocks.
        Returns: (list of 3D concentration grids, per-particle concentrations)
        """
        positions = np.asarray(original_positions)
        mass = np.asarray(mass).ravel()

        if block_indices is None or valid_mask is None:
            x_bounds = np.array([b['x_edges'][0] for b in master_grid] + [master_grid[-1]['x_edges'][-1]])
            block_indices = np.digitize(positions[:, 0], x_bounds) - 1
            block_indices[positions[:, 0] == x_bounds[-1]] = len(master_grid) - 1
            valid_mask = (block_indices >= 0) & (block_indices < len(master_grid))

        results = []
        original_order_concentrations = np.zeros_like(mass)

        for idx, block in enumerate(master_grid):
            x_edges, y_edges, z_edges = block['x_edges'], block['y_edges'], block['z_edges']

            dx = x_edges[1] - x_edges[0]
            dy = y_edges[1] - y_edges[0]
            dz = z_edges[1] - z_edges[0]
            cell_volume = dx * dy * dz

            block_mask = (block_indices == idx) & valid_mask
            if not np.any(block_mask):
                concentration_grid = np.zeros((len(x_edges)-1, len(y_edges)-1, len(z_edges)-1))
            else:
                block_pos = positions[block_mask]
                block_mass = mass[block_mask]

                mass_grid, _ = np.histogramdd(
                    block_pos,
                    bins=[x_edges, y_edges, z_edges],
                    weights=block_mass
                )
                concentration_grid = mass_grid / cell_volume

            results.append(concentration_grid)

            # Lookup for particles in this block
            block_mask = (block_indices == idx) & valid_mask
            if np.any(block_mask):
                block_pos = positions[block_mask]
                ix = np.searchsorted(x_edges, block_pos[:, 0], side='right') - 1
                iy = np.searchsorted(y_edges, block_pos[:, 1], side='right') - 1
                iz = np.searchsorted(z_edges, block_pos[:, 2], side='right') - 1

                ix = np.clip(ix, 0, concentration_grid.shape[0] - 1)
                iy = np.clip(iy, 0, concentration_grid.shape[1] - 1)
                iz = np.clip(iz, 0, concentration_grid.shape[2] - 1)

                original_order_concentrations[block_mask] = concentration_grid[ix, iy, iz]

        return results, original_order_concentrations

    def calculateSecondMoment(self, positions, mass, volume, grid, block_indices=None, valid_mask=None):
        """
        Calculate the second moment of mass across the grid.
        Dispatches to either uniform or non-uniform calculation based on the grid type.

        Args:
            positions (np.ndarray): Array of particle positions [x, y, z].
            mass (np.ndarray): Array of particle masses.
            volume (np.ndarray): Array of particle volumes.
            grid (list|tuple|list[dict]): The grid definition.
            block_indices (np.ndarray, optional): Pre-calculated block assignment.
            valid_mask (np.ndarray, optional): Pre-calculated validity mask.
        """
        is_uniform = (isinstance(grid, (list, tuple)) and len(grid) == 3 and not isinstance(grid[0], dict))
        return self.calculateSecondMoment_uniform(positions, mass, volume, grid) \
            if is_uniform else self.calculateSecondMoment_nonuniform(positions, mass, volume, grid, block_indices, valid_mask)

    def calculateSecondMoment_uniform(self, positions, mass, volume, grid):
        """
        Calculate the second moment for a uniform grid.

        Args:
            positions (np.ndarray): Array of particle positions [x, y, z].
            mass (np.ndarray): Array of particle masses.
            volume (np.ndarray): Array of particle volumes.
            grid (list): A list of 3 arrays containing the edges for x, y, and z.

        Returns:
            np.ndarray: The second moment grid.
        """
        positions = np.asarray(positions)
        mass = np.asarray(mass).ravel()
        volume = np.asarray(volume).ravel()

        # Unpack grid edges
        x_edges, y_edges, z_edges = grid

        # 1. Compute custom weights: mass^2 / volume
        custom_weights = (mass ** 2) / volume

        # 2. 3D Histogram with weights across all dimensions simultaneously
        sum_grid, _ = np.histogramdd(
            positions,
            bins=[x_edges, y_edges, z_edges],
            weights=custom_weights
        )

        # 3. Compute 3D cell volumes for rectilinear grids
        dx = np.diff(x_edges)
        dy = np.diff(y_edges)
        dz = np.diff(z_edges)
        cell_volumes = dx[:, None, None] * dy[None, :, None] * dz[None, None, :]

        # 4. Divide by cell volume to get second moment grid
        second_moment_grid = sum_grid / cell_volumes

        return second_moment_grid

    def calculateSecondMoment_nonuniform(self, positions, mass, volume, master_grid, block_indices=None, valid_mask=None):
        """
        Calculate the second moment for a non-uniform (adaptive) grid using 3D blocks.
        """
        positions = np.asarray(positions)
        mass = np.asarray(mass).ravel()
        volume = np.asarray(volume).ravel()

        if block_indices is None or valid_mask is None:
            x_bounds = np.array([b['x_edges'][0] for b in master_grid] + [master_grid[-1]['x_edges'][-1]])
            block_indices = np.digitize(positions[:, 0], x_bounds) - 1
            block_indices[positions[:, 0] == x_bounds[-1]] = len(master_grid) - 1
            valid_mask = (block_indices >= 0) & (block_indices < len(master_grid))

        second_moment_results = []
        for idx, block in enumerate(master_grid):
            x_edges, y_edges, z_edges = block['x_edges'], block['y_edges'], block['z_edges']
            dx = x_edges[1] - x_edges[0]
            dy = y_edges[1] - y_edges[0]
            dz = z_edges[1] - z_edges[0]
            cell_volume = dx * dy * dz

            block_mask = (block_indices == idx) & valid_mask
            if not np.any(block_mask):
                second_moment_results.append(np.zeros((len(x_edges)-1, len(y_edges)-1, len(z_edges)-1)))
                continue

            block_pos = positions[block_mask]
            # Weights: mass^2 / volume
            weights = (mass[block_mask] ** 2) / volume[block_mask]

            sum_grid, _ = np.histogramdd(
                block_pos,
                bins=[x_edges, y_edges, z_edges],
                weights=weights
            )
            second_moment_results.append(sum_grid / cell_volume)

        return second_moment_results

    def calculateOccupiedVolume(self, positions, volumes, grid, block_indices=None, valid_mask=None):
        """
        Calculate the occupied volume fraction across the grid.
        Dispatches to either uniform or non-uniform calculation based on the grid type.

        Args:
            positions (np.ndarray): Array of particle positions [x, y, z].
            volumes (np.ndarray): Array of particle volumes.
            grid (list|tuple|list[dict]): The grid definition.
            block_indices (np.ndarray, optional): Pre-calculated block assignment.
            valid_mask (np.ndarray, optional): Pre-calculated validity mask.
        """
        is_uniform = (isinstance(grid, (list, tuple)) and len(grid) == 3 and not isinstance(grid[0], dict))
        return self.calculateOccupiedVolume_uniform(positions, volumes, grid) \
            if is_uniform else self.calculateOccupiedVolume_nonuniform(positions, volumes, grid, block_indices, valid_mask)

    def calculateOccupiedVolume_uniform(self, positions, volumes, grid):
        """
        Calculate occupied volume for a uniform grid.

        Args:
            positions (np.ndarray): Array of particle positions [x, y, z].
            volumes (np.ndarray): Array of particle volumes.
            grid (list): A list of 3 arrays containing the edges for x, y, and z.

        Returns:
            np.ndarray: The occupied volume grid.
        """
        x_edges, y_edges, z_edges = grid
        volumes = np.asarray(volumes).ravel()

        # 3D histogram for uniform grid
        sum_vol, _ = np.histogramdd(positions, bins=grid, weights=volumes)

        dx = np.diff(x_edges)
        dy = np.diff(y_edges)
        dz = np.diff(z_edges)
        cell_volumes = dx[:, None, None] * dy[None, :, None] * dz[None, None, :]

        return sum_vol / cell_volumes

    def calculateOccupiedVolume_nonuniform(self, positions, volumes, master_grid, block_indices=None, valid_mask=None):
        """
        Calculate occupied volume for a non-uniform (adaptive) grid using 3D blocks.
        """
        positions = np.asarray(positions)
        volumes = np.asarray(volumes).ravel()

        if block_indices is None or valid_mask is None:
            x_bounds = np.array([b['x_edges'][0] for b in master_grid] + [master_grid[-1]['x_edges'][-1]])
            block_indices = np.digitize(positions[:, 0], x_bounds) - 1
            block_indices[positions[:, 0] == x_bounds[-1]] = len(master_grid) - 1
            valid_mask = (block_indices >= 0) & (block_indices < len(master_grid))

        occupied_volumes_list = []
        for idx, block in enumerate(master_grid):
            x_edges, y_edges, z_edges = block['x_edges'], block['y_edges'], block['z_edges']
            dx = x_edges[1] - x_edges[0]
            dy = y_edges[1] - y_edges[0]
            dz = z_edges[1] - z_edges[0]
            cell_volume = dx * dy * dz

            block_mask = (block_indices == idx) & valid_mask
            if not np.any(block_mask):
                occupied_volumes_list.append(np.zeros((len(x_edges)-1, len(y_edges)-1, len(z_edges)-1)))
                continue

            sum_vol, _ = np.histogramdd(
                positions[block_mask],
                bins=[x_edges, y_edges, z_edges],
                weights=volumes[block_mask]
            )
            occupied_volumes_list.append(sum_vol / cell_volume)

        return occupied_volumes_list

    def adjustTmByCell(self, positions, tm, grid, block_indices=None, valid_mask=None, mode='mean'):
        """
        Adjusts the tm values of particles so that all particles in the same cell have the same tm.
        The value is either the mean or the max of tm values in that cell.

        Args:
            positions (np.ndarray): Array of particle positions [x, y, z].
            tm (np.ndarray): Array of tm values for each particle.
            grid (list|tuple|list[dict]): The grid definition.
            block_indices (np.ndarray, optional): Pre-calculated block assignment.
            valid_mask (np.ndarray, optional): Pre-calculated validity mask.
            mode (str): 'mean' or 'max'.

        Returns:
            np.ndarray: The adjusted tm values.
        """
        if mode == "individual":
            return tm

        is_uniform = (isinstance(grid, (list, tuple)) and len(grid) == 3 and not isinstance(grid[0], dict))
        return self.adjustTmByCell_uniform(positions, tm, grid, mode) \
            if is_uniform else self.adjustTmByCell_nonuniform(positions, tm, grid, block_indices, valid_mask, mode)

    def adjustTmByCell_uniform(self, positions, tm, grid, mode):
        positions = np.asarray(positions)
        tm = np.asarray(tm).ravel()
        x_edges, y_edges, z_edges = grid
        nx, ny, nz = len(x_edges) - 1, len(y_edges) - 1, len(z_edges) - 1

        px, py, pz = positions[:, 0], positions[:, 1], positions[:, 2]
        ix = np.searchsorted(x_edges, px, side='right') - 1
        iy = np.searchsorted(y_edges, py, side='right') - 1
        iz = np.searchsorted(z_edges, pz, side='right') - 1

        valid = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny) & (iz >= 0) & (iz < nz)
        flat_idx = ix[valid] * (ny * nz) + iy[valid] * nz + iz[valid]

        tm_adj = np.copy(tm)
        if not np.any(valid):
            return tm_adj

        if mode == 'mean':
            sums = np.bincount(flat_idx, weights=tm[valid], minlength=nx*ny*nz)
            counts = np.bincount(flat_idx, minlength=nx*ny*nz)
            cell_vals = sums / np.where(counts > 0, counts, 1)
        elif mode == 'max':
            cell_vals = np.full(nx*ny*nz, -np.inf)
            np.maximum.at(cell_vals, flat_idx, tm[valid])
        else:
            raise ValueError("mode must be 'mean' or 'max'")

        tm_adj[valid] = cell_vals[flat_idx]
        return tm_adj

    def adjustTmByCell_nonuniform(self, positions, tm, master_grid, block_indices=None, valid_mask=None, mode='mean'):
        positions = np.asarray(positions)
        tm = np.asarray(tm).ravel()

        if block_indices is None or valid_mask is None:
            x_bounds = np.array([b['x_edges'][0] for b in master_grid] + [master_grid[-1]['x_edges'][-1]])
            block_indices = np.digitize(positions[:, 0], x_bounds) - 1
            block_indices[positions[:, 0] == x_bounds[-1]] = len(master_grid) - 1
            valid_mask = (block_indices >= 0) & (block_indices < len(master_grid))

        tm_adj = np.copy(tm)

        for idx, block in enumerate(master_grid):
            x_edges, y_edges, z_edges = block['x_edges'], block['y_edges'], block['z_edges']
            nx, ny, nz = len(x_edges)-1, len(y_edges)-1, len(z_edges)-1

            block_mask = (block_indices == idx) & valid_mask
            if not np.any(block_mask):
                continue

            b_pos = positions[block_mask]
            b_tm = tm[block_mask]

            ix = np.searchsorted(x_edges, b_pos[:, 0], side='right') - 1
            iy = np.searchsorted(y_edges, b_pos[:, 1], side='right') - 1
            iz = np.searchsorted(z_edges, b_pos[:, 2], side='right') - 1

            b_valid = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny) & (iz >= 0) & (iz < nz)

            if not np.any(b_valid):
                continue

            v_ix = ix[b_valid]
            v_iy = iy[b_valid]
            v_iz = iz[b_valid]
            v_tm = b_tm[b_valid]

            flat_idx = v_ix * (ny * nz) + v_iy * nz + v_iz

            if mode == 'mean':
                sums = np.bincount(flat_idx, weights=v_tm, minlength=nx*ny*nz)
                counts = np.bincount(flat_idx, minlength=nx*ny*nz)
                cell_vals = sums / np.where(counts > 0, counts, 1)
            elif mode == 'max':
                cell_vals = np.full(nx*ny*nz, -np.inf)
                np.maximum.at(cell_vals, flat_idx, v_tm)
            else:
                raise ValueError("mode must be 'mean' or 'max'")

            block_indices_all = np.where(block_mask)[0]
            valid_in_block_indices = block_indices_all[b_valid]

            tm_adj[valid_in_block_indices] = cell_vals[flat_idx]

        return tm_adj

    def updatePhysicsState(self, t0_active, d_r_active, sig_r_active, sigma0_sqrd, Cr, mu, Cphi, dt,
                             particlesEpsilon, age, sigma_sqrd, particlesTL, warn=False):
        """
        Update particle physics state and calculate relaxation time tm.
        """
        # 1. Update t0
        zero_mask = (t0_active == 0)
        t0_active[zero_mask] = (sigma0_sqrd / Cr / particlesEpsilon[zero_mask])**(1. / 3.)

        # 2. Update d_r
        d_r_active += (3 * Cr * particlesEpsilon * (t0_active + age)**2 * dt)

        # 3. Update sigma_r
        new_sigma_r_sqrd = d_r_active / (
            1 + ((d_r_active - sigma0_sqrd) / (sigma0_sqrd + 2 * sigma_sqrd * particlesTL * age))
        )
        sig_r_active = np.maximum(sig_r_active, new_sigma_r_sqrd)

        # 4. Calculate tm
        L_sqrd = (1.5 * sigma_sqrd)**3 / particlesEpsilon**2
        sigma_ur_sqrd = np.minimum(
            sigma_sqrd,
            sigma_sqrd * (sig_r_active / L_sqrd)**(1. / 3.)
        )
        k = 1.5 * sigma_sqrd
        tm = np.minimum(
            mu * np.sqrt(sig_r_active / sigma_ur_sqrd),
            k / particlesEpsilon / Cphi
        )

        if warn and tm.min() < dt:
            print(f"Warning: minimum calculated tm {tm.min()} < dt!")

        tm = np.maximum(tm, dt + 1e-10)

        return t0_active, d_r_active, sig_r_active, tm
