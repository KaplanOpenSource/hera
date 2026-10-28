# Hermes Workflow Toolkit

The `hermesWorkflowToolkit` is the database-backed workflow orchestrator that wraps the [Hermes](https://github.com/KaplanOpenSource/hermes) library. It manages the full lifecycle of simulation workflows: creating them from JSON, storing them in MongoDB with auto-generated names, building them into Luigi task DAGs, executing them, and comparing parameter variations.

`OFToolkit` inherits from `hermesWorkflowToolkit`, gaining workflow support automatically. The base toolkit itself is registered as `toolkitHome.SIMULATIONS_WORKFLOWS` (`"hermesWorkflows"`) and is what the `hera-workflows` CLI uses (handlers in `hera/simulations/CLI.py`).

---

## Architecture overview

### Relationship to Hermes

The toolkit is a **wrapper** around the Hermes workflow engine:

| Layer | Class | Role |
|-------|-------|------|
| **Hera toolkit** | `hermesWorkflowToolkit` | MongoDB storage, naming, retrieval, comparison, execution orchestration |
| **Hermes workflow** | `hermes.workflow` | JSON → task DAG construction, parameter extraction, node management |
| **Hermes engine** | `hermes.engines.luigi.builder` | Task DAG → Luigi Python module code generation |
| **Hermes wrapper** | `hermes.taskwrapper.wrapper` | Wraps each node into a TaskWrapper with dependencies, parameters, and executer resolution |
| **Luigi** | `luigi.Task` | Dependency-based execution engine with target-file state tracking |

### Class hierarchy

```
abstractToolkit
    ↓
hermesWorkflowToolkit (hera/simulations/hermesWorkflowToolkit.py)
    └─ OFToolkit (hera/simulations/openFoam/toolkit.py)
        └─ adds OFObjectHome, Analysis, Presentation, solver extensions
```

`OFToolkit` is the only subclass in master. (`LSMToolkit` inherits from `abstractToolkit` directly; there is no LSM-specific workflow toolkit.)

### Toolkit extension pattern

Solver-specific functionality is added to `OFToolkit` via **toolkit extensions** — composition objects that hold a back-reference to the parent toolkit and provide solver-specific methods. This avoids putting all solver logic into a single class.

```python
# In OFToolkit.__init__:
self.stochasticLagrangian = StochasticLagrangianSolver_toolkitExtension(self)
self.buoyantReactingFoam  = buoyantReactingFoam_toolkitExtension(self)

# Usage:
toolkit.stochasticLagrangian.createDispersionFlowField(...)
toolkit.buoyantReactingFoam.IC_getHydrostaticPressure(...)
```

Each extension follows this pattern:

```python
class absractEulerianSolver_toolkitExtension:
    toolkit = None          # back-reference to parent OFToolkit
    analysis = None         # optional analysis sub-layer
    presentation = None     # optional presentation sub-layer

    def __init__(self, toolkit, solverName, incompressible):
        self.toolkit = toolkit
        self.solverName = solverName
        self.incompressible = incompressible
```

**Extension hierarchy:**

```
absractEulerianSolver_toolkitExtension
    ├─ simpleFoam_toolkitExtension        (not composed into OFToolkit)
    └─ buoyantReactingFoam_toolkitExtension

absractStochasticLagrangianSolver_toolkitExtension
    └─ StochasticLagrangianSolver_toolkitExtension
```

Extensions access the parent toolkit's data layer, workflow management, and mesh utilities via `self.toolkit`. This means they can query MongoDB, read fields, and manage workflows without duplicating any infrastructure:

```python
# Inside an extension method:
docList = self.toolkit.getWorkflowListDocumentFromDB(flowName)
cellCenters = self.toolkit.getMesh(caseDirectory)
self.toolkit.OFObjectHome.readFieldFromCase(...)
```

**Key extension methods:**

| Extension | Method | Purpose |
|-----------|--------|---------|
| `stochasticLagrangian` | `createDispersionFlowField()` | Set up dispersion case from base flow (time mapping, mesh linking) |
| `stochasticLagrangian` | `createDispersionCaseDirectory()` | Validate DB consistency, manage directory conflicts |
| `stochasticLagrangian` | `createAndLinkDispersionCaseDirectory()` | Copy system/constant, symlink processor meshes |
| `stochasticLagrangian` | `writeParticlePositionFile()` | Generate source geometry for particle release |
| `buoyantReactingFoam` | `IC_getHydrostaticPressure()` | Compute hydrostatic pressure field with boundary handling |

### Key constants

```python
DOCTYPE_WORKFLOW = "hermesWorkflow"    # MongoDB document type
DESC_GROUPNAME   = "groupName"         # desc field keys
DESC_GROUPID     = "groupID"
DESC_WORKFLOWNAME = "workflowName"
DESC_PARAMETERS  = "parameters"
DOCKIND_CACHE       = "Cache"          # dockind values for the get*Documents lookups
DOCKIND_SIMULATIONS = "Simulations"

# Module level (hera/simulations/hermesWorkflowToolkit.py):
SCHEDULER_LOCAL   = "local"            # Luigi --local-scheduler
SCHEDULER_CENTRAL = "central"          # connect to a running luigid
class actionModes(Enum): ADD, ADDBUILD, ADDBUILDEXECUTE
```

---

## Workflow JSON structure

A workflow is a directed acyclic graph (DAG) of **nodes**, where each node represents a simulation step (copy files, generate mesh, run solver, etc.):

```json
{
  "workflow": {
    "solver": "simpleFoam",
    "root": "finalnode_xx",
    "nodeList": ["blockMesh", "decomposePar", "simpleFoam", "reconstructPar"],
    "nodes": {
      "blockMesh": {
        "type": "openFOAM.mesh.BlockMesh",
        "Execution": {
          "input_parameters": {
            "vertices": [[0,0,0], [100,0,0], ...],
            "cellCount": [50, 50, 20]
          }
        },
        "requires": []
      },
      "simpleFoam": {
        "type": "openFOAM.system.ControlDict",
        "Execution": {
          "input_parameters": {
            "executeDir": "{blockMesh.output}",
            "endTime": 1000
          }
        },
        "requires": ["blockMesh", "decomposePar"]
      }
    }
  }
}
```

### Parameter referencing

Nodes reference outputs from other nodes using curly-brace syntax:

| Pattern | Meaning |
|---------|---------|
| `{nodeName.param}` | Output parameter from another node |
| `{workflow.param}` | Workflow-level parameter |
| `{WebGui.formData.field}` | GUI form data (Hermes workbench) |

These references are resolved by the `TaskWrapper` to build the dependency graph automatically.

---

## Workflow lifecycle

### 1. Create: load workflow from JSON

```python
# Dynamic class resolution based on solver field.
# If solver is null → hermes.workflow (generic)
# If solver is "simpleFoam" → hera.simulations.openFoam.OFWorkflow.workflow_simpleFoam
wf = toolkit.getHermesWorkflowFromJSON(
    workflow="path/to/workflow.json",  # or dict or JSON string
    name="my_workflow",
    resource="/path/to/case",
)
```

Internally, `pydoc.locate()` dynamically imports the workflow class based on the `solver` field. This allows solver-specific workflow subclasses to add custom node types and validation.

### 2. Add to database

```python
doc = toolkit.addWorkflowToGroup(
    workflowJSON="path/to/workflow.json",
    groupName="dispersion",
    writeWorkflowToFile=False,
    resource=None,                    # auto-generated if None
)
# Creates document: dispersion_0001 (auto-incremented)
```

**Idempotent add:** The method first queries the DB by parameter content. If an identical workflow already exists, it returns the existing document instead of creating a duplicate.

**What gets stored:**
```python
doc = self.addSimulationsDocument(
    resource=resource,
    dataFormat=datatypes.STRING,
    type="hermesWorkflow",
    desc={
        "groupName": groupName,
        "groupID": groupID,
        "workflowName": workflowName,      # e.g. "dispersion_0001"
        "solver": hermesWF.solver,
        "workflow": hermesWF.json,          # full workflow JSON
        "parameters": hermesWF.parametersJSON,  # flattened for querying
    }
)
```

The `parameters` are stored separately from the full `workflow` JSON to enable efficient MongoDB queries via `dictToMongoQuery()` without parsing the workflow tree.

### 3. Build: generate Luigi tasks

```python
# Inside prepareWorkflowRunFromDoc(doc), called by executeWorkflowFromDB():
build = hermesWF.build(buildername=workflow.BUILDER_LUIGI)
```

The build process:

1. **`hermes.workflow._buildNetwork()`** traverses the JSON and creates `TaskWrapper` objects for each node
2. Each `TaskWrapper` extracts its dependencies by scanning `{node.param}` references in input parameters
3. **`LuigiBuilder.buildWorkflow()`** converts the task graph into a Python module containing Luigi Task classes
4. The generated module has one class per task (e.g. `blockMesh_0`, `simpleFoam_0`) with:
   - `requires()` → returns upstream task instances
   - `output()` → returns `luigi.LocalTarget` for state tracking
   - `run()` → calls the Hermes executer for that node type

### 4. Execute: run Luigi pipeline

```python
dispatch_id = toolkit.executeWorkflowFromDB(
    "dispersion_0001",
    scheduler="local",        # or "central" (SCHEDULER_CENTRAL)
    schedulerHost=None, schedulerPort=None,
    dispatch_id=None,         # auto-generated uuid4 hex if None; returned
    stdout=None, stderr=None,
)
```

Execution steps (per matching document; steps 2-5 are `prepareWorkflowRunFromDoc(doc)`):

1. Retrieve workflow document(s) from DB (`getWorkflowListDocumentFromDB`)
2. Reconstruct the workflow object from stored JSON (`getHermesWorkflowFromJSON`)
3. Build Luigi task DAG (generates Python module)
4. Write the workflow JSON to its resource path and the Python module to `filesDirectory/{name}.py`
5. **Clean previous targets** — remove `{name}_targetFiles/` to force re-execution
6. **Run Luigi** via the command from `buildLuigiExecutionCommand()`: `python3 -m luigi --module {name} finalnode_xx_0 --local-scheduler` (local), or with `--scheduler-host/--scheduler-port` for the central scheduler; run with `subprocess.run(..., check=True)`
7. Clean up generated Python module (JSON stays)

The `finalnode_xx` is a synthetic terminal node that depends on all actual nodes — executing it triggers the entire DAG.

### 5. Compare and retrieve

**Cascading search** — `getWorkflowDocumentFromDB(input)` tries four strategies in order:

1. **By name** — exact match on `workflowName` (e.g. `"dispersion_0001"`)
2. **By resource** — match by file/directory path
3. **By group name** — returns all workflows in the group (e.g. `"dispersion"`)
4. **By JSON content** — parses input as JSON, extracts parameters, queries by flattened parameter values via `dictToMongoQuery()`

**Comparison:**

```python
# Compare all workflows in a group
diff = toolkit.compareWorkflowInGroup("dispersion")
# Returns DataFrame: rows = differing parameters, columns = workflow names

# Compare specific workflows
diff = toolkit.compareWorkflows(["dispersion_0001", "dispersion_0002"])

# List all workflows in a group (returns a list of dicts)
toolkit.listWorkflows("dispersion", listNodes=True, listParameters=True)
```

---

## MongoDB storage

### Document schema

```json
{
  "_id": "ObjectId",
  "_cls": "Metadata.Simulations",
  "projectName": "MY_PROJECT",
  "type": "hermesWorkflow",
  "resource": "/path/to/dispersion_0001",
  "dataFormat": "string",
  "desc": {
    "groupName": "dispersion",
    "groupID": 1,
    "workflowName": "dispersion_0001",
    "solver": "simpleFoam",
    "workflow": { "workflow": { "nodeList": [...], "nodes": {...} } },
    "parameters": {
      "blockMesh": { "vertices": [...], "cellCount": [...] },
      "simpleFoam": { "endTime": 1000 }
    }
  }
}
```

### Querying

All queries filter by `type="hermesWorkflow"`. Parameters are queried using `dictToMongoQuery()` which flattens nested dicts into MongoEngine's double-underscore notation:

```python
# {"parameters": {"simpleFoam": {"endTime": 1000}}}
# → {"parameters__simpleFoam__endTime": 1000}
```

---

## Group naming convention

**Format:** `<groupName>_<zero-padded-4-digit-id>`

```python
# Examples: dispersion_0001, flow_0042, lsm_0003
formatted_number = "{0:04d}".format(flowID)
return f"{baseName}_{formatted_number}"
```

**Counter mechanism:** Each group has its own counter stored in the project config via `getCounterAndAdd(groupName)`. The counter auto-increments on each `addWorkflowToGroup()` call that adds a new document, and `deleteWorkflowInGroup(..., resetCounter=True)` resets it. Note that `findAvailableName(simulationGroup)` uses a *different* counter name (`simulations_<group>`), so its result does not track `addWorkflowToGroup`.

**Parsing names:**
```python
base, id = toolkit.splitWorkflowName("dispersion_0001")
# base = "dispersion", id = "0001"
# Underscores inside the base are dropped: "my_group_0001" -> ("mygroup", "0001")
```

---

## CLI (`hera-workflows`)

`hera/bin/hera-workflows` wires argparse subcommands to `hera/simulations/CLI.py`:

| Subcommand | Handler |
|------------|---------|
| `list groups` | `WorkflowsGroup_list` |
| `list workflows` | `workflow_list` |
| `add` | `workflow_add` |
| `buildExecute` | `workflow_buildExecute` |
| `compare` | `workflow_compare` |
| `delete` | `workflow_delete` |
| `export` | `workflow_export` |
| `compareToDisk` | `workflow_compareToDisk` |
| `sync` | `workflow_sync_to_db` |
| `variation` (alias `var`) | `create_workflow_variations` |
| `batch-delete` (alias `bigdelete`) | `batch_delete_workflows` |

See the [CLI reference](../../cli/reference.md#hera-workflows) for arguments.

---

## Key methods reference

### Workflow creation

| Method | Purpose |
|--------|---------|
| `getHermesWorkflowFromJSON(workflow, name=None, resource=None)` | Create workflow object from JSON (file/dict/string); class chosen from `workflow.solver` |
| `getHemresWorkflowFromDocument(documentList, returnFirst)` | Reconstruct workflow from DB document |
| `updateDocumentWorkflow(document, workflow)` | Update stored JSON + parameters |

### Workflow storage

| Method | Purpose |
|--------|---------|
| `addWorkflowToGroup(workflowJSON, groupName, writeWorkflowToFile=False, resource=None)` | Add workflow to DB with auto-naming |
| `addWorkflowFileInGroup(workflowFilePath, write_file=False)` | Add from file (must be inside `FilesDirectory`), infer name/group from filename |
| `findAvailableName(simulationGroup)` | Get next available ID + name |
| `getworkFlowName(baseName, flowID)` | Format `<base>_<padded_id>` |
| `splitWorkflowName(workflow_name)` | Parse `<base>_<id>` |

### Workflow retrieval

| Method | Purpose |
|--------|---------|
| `getHermesWorkflowFromDB(nameOrWorkflowFileOrJSONOrResource, ...)` | Retrieve as hermes.workflow object |
| `getWorkflowDocumentFromDB(nameOrWorkflowFileOrJSONOrResource, doctype=None, dockind="Simulations", **query)` | Retrieve raw document list (cascading search) |
| `getWorkflowDocumentByName(name, doctype=None, dockind="Simulations", **query)` | Retrieve by exact name |
| `getWorkflowListDocumentFromDB(input)` | Retrieve as list of documents |
| `getWorkflowDocumentsInGroup(groupName)` | All documents in a group |
| `getWorkflowListOfSolvers(solverName)` | All documents for a solver |

### Comparison and listing

| Method | Purpose |
|--------|---------|
| `compareWorkflows(Workflow, longFormat=False, transpose=False)` | Compare workflows by name, group, file or list |
| `compareWorkflowInGroup(workflowGroup, longFormat=False, transpose=False)` | Compare all in a group (`None` if the group is empty) |
| `compareWorkflowObj(workflowList, longFormat=False)` | Compare workflow objects directly |
| `listWorkflows(workflowGroup, listNodes=False, listParameters=False)` | List workflows (list of dicts) with optional details |
| `listGroups(solver=None, workflowName=True)` | Print all groups in the project (grouped by solver) |
| `workflowTable(workflowGroup, longFormat=False, transpose=False)` | Alias for `compareWorkflowInGroup` |

### Execution and cleanup

| Method | Purpose |
|--------|---------|
| `executeWorkflowFromDB(nameOrWorkflowFileOrJSONOrResource, scheduler="local", schedulerHost=None, schedulerPort=None, dispatch_id=None, stdout=None, stderr=None)` | Build + execute via Luigi; returns the `dispatch_id` |
| `prepareWorkflowRunFromDoc(doc)` | Build one document's workflow and write the JSON + Luigi module; returns `(pythonFileName, workflowName, targetFilesDir)` |
| `deleteWorkflowInGroup(workflowGroup, deepDelete=False, resetCounter=True, exclude=[])` | Delete workflows, optionally remove files |

### Templates

| Method | Purpose |
|--------|---------|
| `listHermesSolverTemplates(solverName)` | List loaded templates for a solver |
| `getHermesFlowTemplate(hermesFlowName)` | Get a flow template (datasource with `desc.component="Flow"`) |
| `listHermesNodesTemplates()` | List all node templates |
| `getHermesNodeTemplate(hermesNodeName)` | Get a node template (datasource with `desc.component="Node"`) |

---

## Cross-references

| What | Where |
|------|-------|
| User guide | [Toolkits > Simulations > Hermes Workflows](../../toolkits/simulations/workflows.md) |
| OpenFOAM toolkit (inherits from hermesWorkflowToolkit) | [Simulations > OpenFOAM](openfoam.md) |
| LSM toolkit | [Simulations > LSM](lsm.md) |
| Hermes library | `pyHermes/hermes/` (workflow, engines, taskwrapper) |
| API reference | [API > Simulations](../api/simulations.md) |
| Workflow examples | [Examples > Workflows](../../examples/workflows.md) |
