"""Each field the run reads from a sent workflow document.

prepareWorkflowRunFromDoc reads a saved document three ways: doc.desc['workflow'],
doc.desc['workflowName'] and doc['resource']. SentWorkflowDoc wraps the dict the
client sends so those same reads work with no DB lookup. Each test below pins one
of those reads to what was sent.
"""

import pytest

from workflow_runner.sent_workflow_doc import SentWorkflowDoc


def _doc():
    return {
        "desc": {"workflow": {"solver": "s", "nodeList": ["a"]}, "workflowName": "MY_WF"},
        "resource": "/files/MY_WF.json",
    }


def test_desc_workflow_is_the_sent_workflow():
    wrapped = SentWorkflowDoc(_doc())
    # doc.desc['workflow'] -> the workflow JSON the run builds from.
    assert wrapped.desc["workflow"] == {"solver": "s", "nodeList": ["a"]}


def test_desc_workflow_name_is_the_sent_name():
    wrapped = SentWorkflowDoc(_doc())
    # doc.desc['workflowName'] -> the name of the module and target files.
    assert wrapped.desc["workflowName"] == "MY_WF"


def test_resource_is_read_by_item_access():
    wrapped = SentWorkflowDoc(_doc())
    # doc['resource'] -> where the workflow JSON is written (item access, not attribute).
    assert wrapped["resource"] == "/files/MY_WF.json"


def test_matches_every_read_prepare_makes():
    doc = _doc()
    wrapped = SentWorkflowDoc(doc)
    # Mirror the exact three reads prepareWorkflowRunFromDoc makes on a real doc.
    assert wrapped.desc["workflow"] == doc["desc"]["workflow"]
    assert wrapped.desc["workflowName"] == doc["desc"]["workflowName"]
    assert wrapped["resource"] == doc["resource"]


def test_missing_resource_raises_key_error():
    # No silent default: a doc with no resource fails loudly where it is read.
    wrapped = SentWorkflowDoc({"desc": {"workflow": {}, "workflowName": "W"}})
    with pytest.raises(KeyError):
        wrapped["resource"]
