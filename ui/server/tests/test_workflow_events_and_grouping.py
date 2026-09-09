"""Deterministic tests for the two halves of per-task output grouping.

End-to-end grouping through ``luigi.build`` can't be asserted reliably: a single
reader thread tags captured bytes by ``task_pointer.current`` at read time, and a
fast toy DAG changes the pointer faster than the reads, so boundaries race. So we
test the pieces instead:
  1. the real Luigi event handlers move the task pointer, and
  2. WorkflowLogBuilder groups tagged messages into per-task segments.
Together these are the mechanism; the glue between them is a single reader loop.
"""

import types

# Importing this registers the handlers on luigi.Task (needs luigi installed).
import workflow_runner.luigi_task_events  # noqa: F401
from workflow_runner.luigi_task_events import LuigiTaskEvents
from workflow_runner.task_pointer import BETWEEN, PREAMBLE, task_pointer
from workflow_runner.workflow_child_result import WorkflowOutput
from workflow_runner.workflow_log_builder import WorkflowLogBuilder


def _task(family):
    return types.SimpleNamespace(task_family=family)


def test_start_points_at_the_task_then_success_clears_it(capsys):
    task_pointer.current = PREAMBLE
    LuigiTaskEvents.on_start(_task("BuildMesh"))
    assert task_pointer.current == "BuildMesh"
    LuigiTaskEvents.on_success(_task("BuildMesh"))
    assert task_pointer.current == BETWEEN


def test_failure_clears_the_pointer(capsys):
    task_pointer.current = "BuildMesh"
    LuigiTaskEvents.on_failure(_task("BuildMesh"), RuntimeError("boom"))
    assert task_pointer.current == BETWEEN


def test_log_builder_groups_consecutive_output_by_task():
    builder = WorkflowLogBuilder()
    tagged = [
        (PREAMBLE, "scheduling\n"),
        ("TaskA", "a1\n"),
        ("TaskA", "a2\n"),
        ("TaskB", "b\n"),
        (BETWEEN, "done\n"),
    ]
    for task, text in tagged:
        builder.add(WorkflowOutput(task=task, text=text))

    chunks = builder.chunks()

    assert [chunk["name"] for chunk in chunks] == [PREAMBLE, "TaskA", "TaskB", BETWEEN]
    # Consecutive pieces for the same task merge into one segment.
    assert chunks[1]["text"] == "a1\na2\n"
    # The flat log is every piece in arrival order.
    assert builder.output() == "scheduling\na1\na2\nb\ndone\n"


def test_log_builder_keeps_separate_between_segments_in_order():
    builder = WorkflowLogBuilder()
    for task, text in [("TaskA", "a\n"), (BETWEEN, "gap1\n"), ("TaskB", "b\n"), (BETWEEN, "gap2\n")]:
        builder.add(WorkflowOutput(task=task, text=text))

    # Two BETWEEN gaps stay as two segments in run order, not merged into one.
    assert [chunk["name"] for chunk in builder.chunks()] == ["TaskA", BETWEEN, "TaskB", BETWEEN]
