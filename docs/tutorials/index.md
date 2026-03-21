# Tutorials

Interactive Jupyter notebook tutorials that walk you through Hera's core workflows step by step with real code and output.

These notebooks complement the User Guide — they show the same concepts in a hands-on, executable format.

---

## Beginner Tutorials

Start here if you're new to Hera. These notebooks cover the fundamentals in order:

| Notebook | What you'll learn | Related docs | Download |
|----------|------------------|-------------|----------|
| [Project](Project.ipynb) | Creating projects (CLI + Python), caseConfiguration.json, loading documents manually | [Projects](../user_guide/projects.md) | [:material-download: .ipynb](https://github.com/KaplanOpenSource/hera/raw/master/docs/tutorials/Project.ipynb) |
| [Toolkit and DataSource](Toolkit_and_DataSource.ipynb) | What toolkits are, available toolkit names, initializing toolkits, data source JSON structure | [Key Concepts](../user_guide/concepts.md), [Working with Data](../user_guide/working_with_data.md#data-sources) | [:material-download: .ipynb](https://github.com/KaplanOpenSource/hera/raw/master/docs/tutorials/Toolkit_and_DataSource.ipynb) |
| [DataSource](DataSource.ipynb) | Adding data sources, versioning, listing, querying, getting data, deleting | [Working with Data > Data Sources](../user_guide/working_with_data.md#data-sources) | [:material-download: .ipynb](https://github.com/KaplanOpenSource/hera/raw/master/docs/tutorials/DataSource.ipynb) |
| [Repository](Repository.ipynb) | Repository management (add, list, remove, show), JSON structure, DataSource vs Measurements items | [Key Concepts > Repositories](../user_guide/concepts.md#repositories), [CLI Reference](../cli/reference.md#repository-management) | [:material-download: .ipynb](https://github.com/KaplanOpenSource/hera/raw/master/docs/tutorials/Repository.ipynb) |

## Toolkit Tutorials

| Notebook | What you'll learn | Related docs | Download |
|----------|------------------|-------------|----------|
| [Tile Toolkit](TileToolkitDocumentation.ipynb) | Map images from tile servers, WGS84 and ITM coordinates, custom tile servers | [Measurements > Tiles](../toolkits/measurements.md#tiles) | [:material-download: .ipynb](https://github.com/KaplanOpenSource/hera/raw/master/docs/tutorials/TileToolkitDocumentation.ipynb) |

### How to run the tutorials

1. Activate the Hera environment:
   ```bash
   source activate_hera.sh
   ```

2. Make sure MongoDB is running:
   ```bash
   make mongo-up
   ```

3. Start Jupyter:
   ```bash
   jupyter lab docs/tutorials/
   ```

4. Run the notebooks in order: Project → Toolkit and DataSource → DataSource → Repository.

---

## More notebooks

The full collection of Jupyter notebooks (including advanced topics) is in the repository under `hera/doc/jupyter/`:

| Directory | Content |
|-----------|---------|
| `hera/doc/jupyter/User/datalayer/` | Data layer deep dives: low-level access, caching, version management |
| `hera/doc/jupyter/User/toolkits/measurments/GIS/` | GIS toolkits: topography, buildings, demography, land cover, tiles |
| `hera/doc/jupyter/User/toolkits/measurments/experiment/` | Experiment creation and usage |
| `hera/doc/jupyter/User/toolkits/measurments/low_frequency/` | Meteorology low-frequency data |
| `hera/doc/jupyter/User/toolkits/simulations/` | OpenFOAM, Gaussian dispersion, wind profile, Hermes workflows |
| `hera/doc/jupyter/User/toolkits/risskassessment/` | Risk assessment |
| `hera/doc/jupyter/User/utils/` | Utilities: JSON, units (Unum), queries, DataFrames, contours, Slurm |
| `hera/doc/jupyter/Developer/` | Developer-focused notebooks with additional low-level details |
