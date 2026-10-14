"""Objects the routers share, created once for the process."""
from hera_warmup import HeraWarmup
from workflow_runner import WorkflowRunner

# True once hera finished its slow imports; every route that reaches hera waits for it.
warmup = HeraWarmup()

# Holds the state of the one workflow run the server allows at a time.
workflow_runner = WorkflowRunner()
