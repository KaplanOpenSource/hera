import numpy as np

class volumetricParticleCalculator:


    def calculateConcentration(self, original_positions, mass, grid):
        """
        Calculate the mass concentration across the grid.
        Dispatches to either uniform or non-uniform calculation based on the grid type.

        Args:
            original_positions (np.ndarray): Array of particle positions [x, y, z].
            mass (np.ndarray): Array of particle masses.
            grid (list|tuple|list[dict]): The grid definition (edges for uniform, master_grid for non-uniform).

        Returns:
            tuple: (concentration_grid, original_order_concentrations)
        """
        is_uniform = (isinstance(grid, (list, tuple)) and len(grid) == 3 and not isinstance(grid[0], dict))
        return self.calculateConcentration_uniform(original_positions, mass, grid) \
            if is_uniform else self.calculateConcentration_nonuniform(original_positions, mass, grid)

    def calculateConcentration_uniform(self, original_positions, mass, grid):
        """
        Calculate concentration for a uniform grid using histogramdd.

        Args:
            original_positions (np.ndarray): Array of particle positions [x, y, z].
            mass (np.ndarray): Array of particle masses.
            grid (list): A list of 3 arrays containing the edges for x, y, and z.

        Returns:
            tuple: (concentration_grid, original_order_concentrations)
        """
        original_positions = np.asarray(original_positions)
        mass = np.asarray(mass)
        
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

    def calculateConcentration_nonuniform(self, original_positions, mass, master_grid):
        """
        Calculate concentration for a non-uniform (adaptive) grid using 3D blocks.
        Returns: (list of 3D concentration grids, per-particle concentrations)
        """
        positions = np.asarray(original_positions)
        mass = np.asarray(mass)

        results = []
        original_order_concentrations = np.zeros_like(mass)

        for idx, block in enumerate(master_grid):
            x_edges, y_edges, z_edges = block['x_edges'], block['y_edges'], block['z_edges']
            dx = x_edges[1] - x_edges[0]
            dy = y_edges[1] - y_edges[0]
            dz = z_edges[1] - z_edges[0]
            cell_volume = dx * dy * dz

            # 1. Filter particles in this block
            # Using the same boundaries as get_adaptive_grid (inclusive end for last block)
            if idx == len(master_grid) - 1:
                mask = (positions[:, 0] >= x_edges[0]) & (positions[:, 0] <= x_edges[-1])
            else:
                mask = (positions[:, 0] >= x_edges[0]) & (positions[:, 0] < x_edges[-1])

            block_pos = positions[mask]
            block_mass = mass[mask]

            # 2. 3D Histogram for the block
            mass_grid, _ = np.histogramdd(
                block_pos,
                bins=[x_edges, y_edges, z_edges],
                weights=block_mass
            )

            concentration_grid = mass_grid / cell_volume
            results.append(concentration_grid)

            # 3. Look up concentrations for particles in this block
            if len(block_pos) > 0:
                # Vectorized index lookup
                ix = np.searchsorted(x_edges, block_pos[:, 0], side='right') - 1
                iy = np.searchsorted(y_edges, block_pos[:, 1], side='right') - 1
                iz = np.searchsorted(z_edges, block_pos[:, 2], side='right') - 1

                # Clamp to avoid boundary errors
                ix = np.clip(ix, 0, concentration_grid.shape[0] - 1)
                iy = np.clip(iy, 0, concentration_grid.shape[1] - 1)
                iz = np.clip(iz, 0, concentration_grid.shape[2] - 1)

                original_order_concentrations[mask] = concentration_grid[ix, iy, iz]

        return results, original_order_concentrations

    def calculateSecondMoment(self, positions, mass, volume, grid):
        """
        Calculate the second moment of mass across the grid.
        Dispatches to either uniform or non-uniform calculation based on the grid type.

        Args:
            positions (np.ndarray): Array of particle positions [x, y, z].
            mass (np.ndarray): Array of particle masses.
            volume (np.ndarray): Array of particle volumes.
            grid (list|tuple|list[dict]): The grid definition.

        Returns:
            np.ndarray|list: The second moment grid (or list of grids for non-uniform).
        """
        is_uniform = (isinstance(grid, (list, tuple)) and len(grid) == 3 and not isinstance(grid[0], dict))
        return self.calculateSecondMoment_uniform(positions, mass, volume, grid) \
            if is_uniform else self.calculateSecondMoment_nonuniform(positions, mass, volume, grid)

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
        mass = np.asarray(mass)
        volume = np.asarray(volume)
        
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

    def calculateSecondMoment_nonuniform(self, positions, mass, volume, master_grid):
        """
        Calculate the second moment for a non-uniform (adaptive) grid using 3D blocks.
        """
        positions = np.asarray(positions)
        mass = np.asarray(mass)
        volume = np.asarray(volume)

        second_moment_results = []
        for idx, block in enumerate(master_grid):
            x_edges, y_edges, z_edges = block['x_edges'], block['y_edges'], block['z_edges']
            dx = x_edges[1] - x_edges[0]
            dy = y_edges[1] - y_edges[0]
            dz = z_edges[1] - z_edges[0]
            cell_volume = dx * dy * dz

            if idx == len(master_grid) - 1:
                mask = (positions[:, 0] >= x_edges[0]) & (positions[:, 0] <= x_edges[-1])
            else:
                mask = (positions[:, 0] >= x_edges[0]) & (positions[:, 0] < x_edges[-1])

            block_pos = positions[mask]
            # Weights: mass^2 / volume
            weights = (mass[mask] ** 2) / volume[mask]

            sum_grid, _ = np.histogramdd(
                block_pos,
                bins=[x_edges, y_edges, z_edges],
                weights=weights
            )
            second_moment_results.append(sum_grid / cell_volume)

        return second_moment_results

    def calculateOccupiedVolume(self, positions, volumes, grid):
        """
        Calculate the occupied volume fraction across the grid.
        Dispatches to either uniform or non-uniform calculation based on the grid type.

        Args:
            positions (np.ndarray): Array of particle positions [x, y, z].
            volumes (np.ndarray): Array of particle volumes.
            grid (list|tuple|list[dict]): The grid definition.

        Returns:
            np.ndarray|list: The occupied volume grid (or list of grids for non-uniform).
        """
        is_uniform = (isinstance(grid, (list, tuple)) and len(grid) == 3 and not isinstance(grid[0], dict))
        return self.calculateOccupiedVolume_uniform(positions, volumes, grid) \
            if is_uniform else self.calculateOccupiedVolume_nonuniform(positions, volumes, grid)

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
        
        # 3D histogram for uniform grid
        sum_vol, _ = np.histogramdd(positions, bins=grid, weights=volumes)
        
        dx = np.diff(x_edges)
        dy = np.diff(y_edges)
        dz = np.diff(z_edges)
        cell_volumes = dx[:, None, None] * dy[None, :, None] * dz[None, None, :]
        
        return sum_vol / cell_volumes

    def calculateOccupiedVolume_nonuniform(self, positions, volumes, master_grid):
        """
        Calculate occupied volume for a non-uniform (adaptive) grid using 3D blocks.
        """
        positions = np.asarray(positions)
        volumes = np.asarray(volumes)

        occupied_volumes_list = []
        for idx, block in enumerate(master_grid):
            x_edges, y_edges, z_edges = block['x_edges'], block['y_edges'], block['z_edges']
            dx = x_edges[1] - x_edges[0]
            dy = y_edges[1] - y_edges[0]
            dz = z_edges[1] - z_edges[0]
            cell_volume = dx * dy * dz

            if idx == len(master_grid) - 1:
                mask = (positions[:, 0] >= x_edges[0]) & (positions[:, 0] <= x_edges[-1])
            else:
                mask = (positions[:, 0] >= x_edges[0]) & (positions[:, 0] < x_edges[-1])

            sum_vol, _ = np.histogramdd(
                positions[mask],
                bins=[x_edges, y_edges, z_edges],
                weights=volumes[mask]
            )
            occupied_volumes_list.append(sum_vol / cell_volume)

        return occupied_volumes_list
