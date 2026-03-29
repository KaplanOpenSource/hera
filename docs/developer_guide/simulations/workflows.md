# Hermes Workflow Toolkit

The `hermesWorkflowToolkit` is the database-backed workflow orchestrator that wraps the [Hermes](https://github.com/KaplanOpenSource/hermes) library. It manages the full lifecycle of simulation workflows: creating them from JSON, storing them in MongoDB with auto-generated names, building them into Luigi task DAGs, executing them, and comparing parameter variations.

`OFToolkit` inherits from `hermesWorkflowToolkit`, gaining workflow support automatically.

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
    ├─ OFToolkit (hera/simulations/openFoam/toolkit.py)
    │   └─ adds OFObjectHome, Analysis, Presentation, solver extensions
    └─ workflowToolkit [LSM variant] (hera/simulations/LSM/hermesWorkflowToolkit.py)
        └─ imports hermes handlers (handler_build, handler_expand, handler_execute)
```

### Key constants

```python
DOCTYPE_WORKFLOW = "hermesWorkflow"    # MongoDB document type
DESC_GROUPNAME   = "groupName"         # desc field keys
DESC_GROUPID     = "groupID"
DESC_WORKFLOWNAME = "workflowName"
DESC_PARAMETERS  = "parameters"
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
# Inside executeWorkflowFromDB():
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
toolkit.executeWorkflowFromDB("dispersion_0001")
```

Execution steps:

1. Retrieve workflow document from DB
2. Reconstruct `hermes.workflow` object from stored JSON
3. Build Luigi task DAG (generates Python module)
4. Write workflow JSON + Python module to `filesDirectory`
5. **Clean previous targets** — remove `{name}_targetFiles/` to force re-execution
6. **Run Luigi** via command line: `python3 -m luigi --module {name} finalnode_xx_0 --local-scheduler`
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

# List all workflows in a group
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

**Counter mechanism:** Each group has its own counter stored in the project config via `getCounterAndAdd(groupName)`. The counter auto-increments on each `addWorkflowToGroup()` call.

**Parsing names:**
```python
base, id = toolkit.splitWorkflowName("dispersion_0001")
# base = "dispersion", id = "0001"
```

---

## LSM variant

**Source:** `hera/simulations/LSM/hermesWorkflowToolkit.py`

The LSM `workflowToolkit` extends the base with additional handler imports for the expand → build → execute pipeline:

| Handler | Purpose |
|---------|---------|
| `handler_expand` | Expand workflow templates with meteorological data |
| `handler_build` | Build the workflow into executable tasks |
| `handler_buildExecute` | Combined build + execute in one step |
| `handler_execute` | Execute a pre-built workflow |

This supports LSM's pattern where workflows are first **expanded** with site-specific meteorological data before being built and executed. The base toolkit only uses build → execute.

---

## Key methods reference

### Workflow creation

| Method | Purpose |
|--------|---------|
| `getHermesWorkflowFromJSON(workflow, name, resource)` | Create workflow object from JSON (file/dict/string) |
| `getHemresWorkflowFromDocument(documentList, returnFirst)` | Reconstruct workflow from DB document |
| `updateDocumentWorkflow(document, workflow)` | Update stored JSON + parameters |

### Workflow storage

| Method | Purpose |
|--------|---------|
| `addWorkflowToGroup(workflowJSON, groupName, ...)` | Add workflow to DB with auto-naming |
| `addWorkflowFileInGroup(workflowFilePath, write_file)` | Add from file, infer group from filename |
| `findAvailableName(simulationGroup)` | Get next available ID + name |
| `getworkFlowName(baseName, flowID)` | Format `<base>_<padded_id>` |
| `splitWorkflowName(workflow_name)` | Parse `<base>_<id>` |

### Workflow retrieval

| Method | Purpose |
|--------|---------|
| `getHermesWorkflowFromDB(input, returnFirst)` | Retrieve as hermes.workflow object |
| `getWorkflowDocumentFromDB(input, doctype, dockind)` | Retrieve raw document (cascading search) |
| `getWorkflowDocumentByName(name, doctype, dockind)` | Retrieve by exact name |
| `getWorkflowListDocumentFromDB(input)` | Retrieve as list of documents |
| `getWorkflowDocumentsInGroup(groupName)` | All documents in a group |
| `getWorkflowListOfSolvers(solverName)` | All documents for a solver |

### Comparison and listing

| Method | Purpose |
|--------|---------|
| `compareWorkflows(Workflow, longFormat, transpose)` | Compare workflows by name or list |
| `compareWorkflowInGroup(workflowGroup)` | Compare all in a group |
| `compareWorkflowObj(workflowList)` | Compare workflow objects directly |
| `listWorkflows(workflowGroup, listNodes, listParameters)` | List workflows with optional details |
| `listGroups(solver, workflowName)` | List all groups in the project |
| `workflowTable(workflowGroup)` | Alias for `compareWorkflowInGroup` |

### Execution and cleanup

| Method | Purpose |
|--------|---------|
| `executeWorkflowFromDB(input)` | Build + execute via Luigi |
| `deleteWorkflowInGroup(workflowGroup, deepDelete, resetCounter)` | Delete workflows, optionally remove files |

### Templates

| Method | Purpose |
|--------|---------|
| `listHermesSolverTemplates(solverName)` | List loaded templates for a solver |
| `getHermesFlowTemplate(hermesFlowName)` | Get a flow template |
| `listHermesNodesTemplates()` | List all node templates |
| `getHermesNodeTemplate(hermesNodeName)` | Get a node template |

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
