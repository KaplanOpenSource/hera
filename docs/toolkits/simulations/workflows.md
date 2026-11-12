# Hermes Workflows

**Toolkit name:** `hermesWorkflows` | **Constant:** `toolkitHome.SIMULATIONS_WORKFLOWS`

The Hermes workflow toolkit stores [Hermes](https://github.com/KaplanOpenSource/hermes) workflows (JSON pipelines of nodes) in the project's **Simulations** collection, organises them into **groups**, compares their parameters, and builds and executes them as [Luigi](https://luigi.readthedocs.io) task DAGs. It is the base class of the [OpenFOAM toolkit](openfoam.md), so everything on this page also works on the OpenFOAM toolkit object.

```python
from hera import toolkitHome

wftk = toolkitHome.getToolkit(toolkitHome.SIMULATIONS_WORKFLOWS, projectName="MY_PROJECT")
```

**Requirements:** the `hermes` package (needed to add, compare, build and execute workflows; adding and comparing raise `NotImplementedError` without it) and `luigi` for execution. Workflow JSON files and the generated Luigi modules are written to the toolkit's files directory (`filesDirectory` argument of `getToolkit`; otherwise the `filesDirectory` saved in the project config, or `~/.hera/<projectName>`).

---

## Groups and names

Every stored workflow belongs to a group and is named `<groupName>_<NNNN>`, where `NNNN` is a 4-digit counter kept per group in the project (e.g. `wind_study_0001`, `wind_study_0002`). The document is of type `hermesWorkflow` and its description holds:

| Key | Content |
|-----|---------|
| `groupName`, `groupID` | The group and the counter value |
| `workflowName` | `<groupName>_<NNNN>` |
| `solver` | The `solver` field of the workflow |
| `workflow` | The full workflow JSON |
| `parameters` | The node parameters, stored separately for querying |

The document resource is the path of the workflow (by default `<filesDirectory>/<workflowName>.json`).

Adding is **idempotent**: if a workflow with identical parameters is already in the DB, the existing document is returned instead of a new one.

---

## Adding workflows

```python
# Add a workflow JSON (file name, JSON string or dict) to a group
doc = wftk.addWorkflowToGroup(
    workflowJSON="/path/to/workflow.json",
    groupName="wind_study",
    writeWorkflowToFile=True,   # also write <filesDirectory>/wind_study_0001.json
)
print(doc.desc["workflowName"])   # wind_study_0001

# Add a file whose name follows <group>_<id>.json; the group is taken from the file name.
# The file must be inside the toolkit files directory.
doc = wftk.addWorkflowFileInGroup("wind_study_1.json")
```

---

## Retrieving workflows

Most retrieval methods accept a workflow name, a group name, a workflow file, a workflow dict, or the resource (case directory):

```python
docs = wftk.getWorkflowListDocumentFromDB("wind_study_0001")   # list of documents (a list of names is also accepted)
doc = wftk.getWorkflowDocumentByName("wind_study_0001")       # by exact name, or None
docs = wftk.getWorkflowDocumentsInGroup("wind_study")          # all documents of a group
docs = wftk.getWorkflowListOfSolvers("simpleFoam")             # all documents of a solver

wf = wftk.getHermesWorkflowFromDB("wind_study_0001")           # hermes workflow object
wf = wftk.getHermesWorkflowFromJSON("/path/to/workflow.json")  # from a file/JSON/dict

# After changing a workflow object, store it back
wftk.updateDocumentWorkflow(document=doc, workflow=wf)
```

The class of the returned workflow object is chosen from its `solver` field (for OpenFOAM solvers, e.g. `hera.simulations.openFoam.OFWorkflow.workflow_simpleFoam`).

---

## Listing and comparing

```python
# Print the groups of the project (per solver), optionally with their workflow names
wftk.listGroups(solver=None, workflowName=True)

# List the workflows of a group: list of dicts with workflowName (+ nodes / parameters)
wftk.listWorkflows("wind_study", listNodes=True, listParameters=False)

# DataFrame of the parameters that differ between workflows
wftk.compareWorkflows(["wind_study_0001", "wind_study_0002"], longFormat=False, transpose=False)
wftk.compareWorkflowInGroup("wind_study")     # all workflows of a group (same as workflowTable)
```

---

## Executing

```python
dispatch_id = wftk.executeWorkflowFromDB(
    "wind_study_0001",
    scheduler="local",        # "local" (luigi --local-scheduler) or "central" (a running luigid)
    schedulerHost=None,       # central scheduler address (default localhost)
    schedulerPort=None,       # (default 8082)
    dispatch_id=None,         # generated (uuid4) if not given; returned
)
```

For every matching document this rebuilds the workflow from the stored JSON, writes the workflow JSON and a generated Python module `<filesDirectory>/<workflowName>.py`, removes `<workflowName>_targetFiles` (so all tasks rerun), runs

```bash
python3 -m luigi --module <workflowName> finalnode_xx_0 --local-scheduler
```

in a subprocess and then deletes the generated module. `prepareWorkflowRunFromDoc(doc)` performs only the preparation steps.

For OpenFOAM workflows use `OFToolkit.runOFSimulation`, which also runs the case's `Allrun` script.

---

## Deleting

```python
# Delete all workflows of a group; deepDelete also removes their resources from disk
wftk.deleteWorkflowInGroup("wind_study", deepDelete=False, resetCounter=True, exclude=["wind_study_0001"])
```

---

## Templates

Workflow and node templates are data sources of the toolkit (loaded from a repository):

```python
wftk.listHermesSolverTemplates("simpleFoam")   # DataFrame of the templates of a solver
wftk.getHermesFlowTemplate("myFlow")            # data source with desc component "Flow"
wftk.listHermesNodesTemplates()                 # node templates (desc component "Node")
wftk.getHermesNodeTemplate("blockMesh")
```

---

## CLI

`hera-workflows` operates on the workflows of a project. When `--projectName` is omitted it is read from `caseConfiguration.json` in the current directory.

```bash
# List the groups (per solver); --workflowName also lists the workflows of each group
hera-workflows list groups --projectName MY_PROJECT [--solver simpleFoam] [--workflowName]

# List the workflows in a group (the argument can also be a workflow name of the group)
hera-workflows list workflows wind_study --projectName MY_PROJECT [--nodes] [--parameters]

# Add a workflow file (named <group>_<id>.json, inside the files directory) to the DB
hera-workflows add wind_study_1.json --projectName MY_PROJECT

# Add (if needed) and execute a workflow; the argument is a DB name or a workflow file
hera-workflows buildExecute wind_study_0001 --projectName MY_PROJECT \
    [--scheduler local|central] [--scheduler-host HOST] [--scheduler-port PORT] [--dispatch-id ID]

# Compare workflows (names, group names or files)
hera-workflows compare wind_study --projectName MY_PROJECT \
    [--longTable] [--transpose] [--format pandas|json|latex|csv] [--file out]

# Export workflows to <workflowName>.json in the current directory
hera-workflows export wind_study_0001 --projectName MY_PROJECT [--overwrite]

# Show the parameter differences between the DB and the file on disk
hera-workflows compareToDisk wind_study_0001 --projectName MY_PROJECT

# Update the DB from the file on disk when the parameters changed (--force: always)
hera-workflows sync wind_study_0001 --projectName MY_PROJECT [--force] [--quiet]

# Delete workflows from the DB (--export writes them to disk first)
hera-workflows delete wind_study_0001 --projectName MY_PROJECT [--export] [--overwrite]

# Create variations of a workflow from a variation JSON and add them to its group
hera-workflows variation wind_study_0001 variations.json --projectName MY_PROJECT [--dry-run] [--overwrite]

# Delete all workflows of a group (except --workflow); --hard also deletes their directories
hera-workflows batch-delete --groupName wind_study --projectName MY_PROJECT [--workflow wind_study_0001] [--hard]
```

Notes:

- `variation` (alias `var`) writes `<group>_<i>.json` next to the base workflow, for every combination of the variation file (a list of variation groups, see `hera.utils.jsonutils.JSONVariations`), and adds each one with the same logic as `add`.
- `delete` writes a `completeDelete.py` script with the commands that remove the deleted workflows' directories (the message it prints calls it `completeRemove.py`).
- `add` accepts `--workflowGroup`, `--overwrite`, `--assignName`, `--allowDuplicate` and `--execute`, but in master they are ignored: the group is always inferred from the file name.

For implementation details, see the [Developer Guide](../../developer_guide/simulations/workflows.md). For the full API, see the [API Reference](../../developer_guide/api/simulations.md).
