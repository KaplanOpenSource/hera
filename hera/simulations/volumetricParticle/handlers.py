import h5py

class DataHandler_adaptive_grid(object):           
    """                                             
    Handler for the adaptive grid structure used in volumetric particle analysis.           
    The resource is a path to an HDF5 file.           
    """                     

    @staticmethod
    def saveData(resource, fileName, **kwargs):
        """
        Save the adaptive grid structure to HDF5.
        'resource' is the list of grid slices (the data to be saved).
        """
        with h5py.File(fileName, "w") as f:
            for idx, slice_info in enumerate(resource):
                grp = f.create_group(f"slice_{idx:03d}")

                # Support both new 'x_edges' and old 'x_range' formats
                x_edges = slice_info.get('x_edges')
                x_range = slice_info.get('x_range')

                if x_edges is not None:
                    grp.create_dataset('x_edges', data=x_edges, compression="gzip")
                    grp.attrs['x_start'] = x_edges[0]
                    grp.attrs['x_end'] = x_edges[-1]
                elif x_range is not None:
                    grp.attrs['x_start'] = x_range[0]
                    grp.attrs['x_end'] = x_range[1]

                grp.create_dataset('y_edges', data=slice_info['y_edges'], compression="gzip")
                grp.create_dataset('z_edges', data=slice_info['z_edges'], compression="gzip")
        return dict()

    @staticmethod
    def getData(resource, desc=None, **kwargs):
        """
        Load the adaptive grid structure from HDF5.
        'resource' is the path to the HDF5 file.
        """

        master_grid = []
        with h5py.File(resource, "r") as f:
            for slice_name in sorted(f.keys()):
                grp = f[slice_name]

                # Load x_edges if available, otherwise fallback to x_range
                if 'x_edges' in grp:
                    x_edges = grp['x_edges'][:]
                    x_range = (x_edges[0], x_edges[-1])
                else:
                    x_edges = None
                    x_range = (grp.attrs['x_start'], grp.attrs['x_end'])

                master_grid.append({
                    'x_range': x_range,
                    'x_edges': x_edges,
                    'y_edges': grp['y_edges'][:],
                    'z_edges': grp['z_edges'][:]
                })
        return sorted(master_grid, key=lambda s: s['x_range'][0]) 

class DataHandler_uniform_grid(object):
    """
    Handler for uniform grid structures.
    The resource is a path to an HDF5 file containing x, y, and z edges.
    """
    @staticmethod
    def saveData(resource, fileName, **kwargs):
        
        # resource is [x_edges, y_edges, z_edges]
        with h5py.File(fileName, "w") as f:
            f.create_dataset('x_edges', data=resource[0], compression="gzip")
            f.create_dataset('y_edges', data=resource[1], compression="gzip")
            f.create_dataset('z_edges', data=resource[2], compression="gzip")
        return dict()

    @staticmethod
    def getData(resource, desc=None, **kwargs):
        
        with h5py.File(resource, "r") as f:
            return [f['x_edges'][:], f['y_edges'][:], f['z_edges'][:]]

class DataHandler_max_cell_grids(object):
    """
    Handler for the maximum cell grids computed over time.
    The resource is a path to an HDF5 file.
    """
    @staticmethod
    def saveData(resource, fileName, **kwargs):
        # resource is a list of dicts: [{'concentration': arr, 'second_moment': arr}, ...]
        with h5py.File(fileName, "w") as f:
            for idx, grid_dict in enumerate(resource):
                grp = f.create_group(f"slice_{idx:03d}")
                if 'concentration' in grid_dict:
                    grp.create_dataset('concentration', data=grid_dict['concentration'], compression="gzip")
                if 'second_moment' in grid_dict:
                    grp.create_dataset('second_moment', data=grid_dict['second_moment'], compression="gzip")
                if 'occupied_volume' in grid_dict:
                    grp.create_dataset('occupied_volume', data=grid_dict['occupied_volume'], compression="gzip")
                if 'fluctuation_intensity' in grid_dict:
                    grp.create_dataset('fluctuation_intensity', data=grid_dict['fluctuation_intensity'], compression="gzip")
                if 'coefficient_of_variation' in grid_dict:
                    grp.create_dataset('coefficient_of_variation', data=grid_dict['coefficient_of_variation'], compression="gzip")
        return dict()

    @staticmethod
    def getData(resource, desc=None, **kwargs):

        max_results = []
        with h5py.File(resource, "r") as f:
            for slice_name in sorted(f.keys()):
                grp = f[slice_name]
                slice_data = {}
                if 'concentration' in grp:
                    slice_data['concentration'] = grp['concentration'][:]
                if 'second_moment' in grp:
                    slice_data['second_moment'] = grp['second_moment'][:]
                if 'occupied_volume' in grp:
                    slice_data['occupied_volume'] = grp['occupied_volume'][:]
                if 'fluctuation_intensity' in grp:
                    slice_data['fluctuation_intensity'] = grp['fluctuation_intensity'][:]
                if 'coefficient_of_variation' in grp:
                    slice_data['coefficient_of_variation'] = grp['coefficient_of_variation'][:]
                max_results.append(slice_data)
        return max_results