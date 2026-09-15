import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { act, cleanup, fireEvent, render, screen } from '@testing-library/react';

// The editor reads the node catalog from a zustand store and mounts a reader
// that fetches it; stub both so nothing hits the network.
vi.mock('../src/components/workflow/useNodeCatalog', () => ({
  useNodeCatalog: (selector: (state: { catalog: unknown[] }) => unknown) => selector({ catalog: [] }),
  NodeCatalogReader: () => null,
}));

// Replace the ReactFlow-based graph with a stub that captures its props, so a
// node edit can be driven without rendering ReactFlow. The autofill still runs
// through the real chain: the editor -> DetailsViewDocument -> fillProjectName.
let graphProps: any = null;
vi.mock('../src/components/workflow/WorkflowGraph', () => ({
  WorkflowGraph: (props: any) => {
    graphProps = props;
    return null;
  },
}));

import { DetailsViewDocument } from '../src/components/details/DetailsViewDocument';
import { DocumentObj, ProjectObj } from '../src/objects/ProjectObj';
import { fillProjectName, isProjectNameKey } from '../src/shared/workflowMutators/mutators/FillProjectNameMutator';
import { WORKFLOW_DOC_TYPE } from '../src/shared/workflow';
import { NO_PROJECT, useProjectStore } from '../src/stores/useProjectStore';

const PROJECT = 'MY_PROJECT';

describe('isProjectNameKey', () => {
  it('matches ProjectName in any casing', () => {
    expect(isProjectNameKey('ProjectName')).toBe(true);
    expect(isProjectNameKey('projectName')).toBe(true);
    expect(isProjectNameKey('projectname')).toBe(true);
    expect(isProjectNameKey('PROJECTNAME')).toBe(true);
  });

  it('does not match other keys', () => {
    expect(isProjectNameKey('project')).toBe(false);
    expect(isProjectNameKey('projectNames')).toBe(false);
    expect(isProjectNameKey('SimulationName')).toBe(false);
  });
});

describe('fillProjectName', () => {
  it('fills a field directly under desc', () => {
    expect(fillProjectName({ projectname: '' }, PROJECT)).toEqual({ projectname: PROJECT });
  });

  it('fills a field nested deep inside node parameters', () => {
    const desc = { workflow: { nodes: { A: { Execution: { input_parameters: { ProjectName: '' } } } } } };
    const result = fillProjectName(desc, PROJECT);
    expect(result.workflow.nodes.A.Execution.input_parameters.ProjectName).toBe(PROJECT);
  });

  it('fills inside arrays', () => {
    const desc = { steps: [{ ProjectName: '' }, { ProjectName: 'KEEP' }] };
    expect(fillProjectName(desc, PROJECT).steps).toEqual([{ ProjectName: PROJECT }, { ProjectName: 'KEEP' }]);
  });

  it('fills a null or unset value', () => {
    const result = fillProjectName({ a: { ProjectName: null }, b: { ProjectName: undefined } }, PROJECT);
    expect(result.a.ProjectName).toBe(PROJECT);
    expect(result.b.ProjectName).toBe(PROJECT);
  });

  it('keeps a value the user set to another project', () => {
    expect(fillProjectName({ ProjectName: 'OTHER' }, PROJECT)).toEqual({ ProjectName: 'OTHER' });
  });

  it('never adds the field', () => {
    expect(fillProjectName({ Command: 'ls' }, PROJECT)).toEqual({ Command: 'ls' });
  });

  it('leaves other values as they are', () => {
    const desc = { n: 0, b: false, s: 'x', nil: null, list: [1, 2] };
    expect(fillProjectName(desc, PROJECT)).toEqual(desc);
  });

  it('does not mutate the input', () => {
    const desc = { ProjectName: '' };
    fillProjectName(desc, PROJECT);
    expect(desc).toEqual({ ProjectName: '' });
  });
});

// A workflow document holding one node, rendered by the workflow editor.
const makeWorkflowDoc = (params: { [key: string]: any }, descExtra: { [key: string]: any } = {}) => {
  const data = {
    _id: { $oid: '1' },
    _cls: 'C',
    type: WORKFLOW_DOC_TYPE,
    resource: '/r',
    desc: {
      ...descExtra,
      workflow: {
        nodeList: ['A'],
        nodes: { A: { type: 'general.RunOsCommand', Execution: { input_parameters: params } } },
      },
    },
  } as any;
  const project = new ProjectObj({ name: PROJECT, documents: [data] } as any);
  return new DocumentObj(data, project);
};

const shownParams = () => {
  return graphProps.nodes.A.Execution.input_parameters;
};

// Writes node A back through the editor, the way the graph does on an edit.
const editNode = (params: { [key: string]: any }) => {
  act(() => graphProps.onSetNode('A', {
    ...graphProps.nodes.A,
    Execution: { ...graphProps.nodes.A.Execution, input_parameters: params },
  }));
};

// Renames a field row in the details tree, the way a user does: click the name
// to start editing, type, then blur to save.
const renameField = (from: string, to: string) => {
  fireEvent.click(screen.getByText(from));
  const input = screen.getByDisplayValue(from);
  fireEvent.change(input, { target: { value: to } });
  fireEvent.blur(input);
};

beforeEach(() => {
  graphProps = null;
  useProjectStore.getState().selectProject(PROJECT);
});

afterEach(() => cleanup());

describe('ProjectName autofill in the document view', () => {
  it('fills a desc field renamed to projectname', () => {
    render(<DetailsViewDocument doc={makeWorkflowDoc({ Command: 'ls' }, { newItem_1: '' })} setDoc={vi.fn()} />);

    renameField('newItem_1', 'projectname');

    expect(screen.getByDisplayValue(PROJECT)).toBeDefined();
  });

  it('fills a ProjectName the user just added to a node', () => {
    render(<DetailsViewDocument doc={makeWorkflowDoc({ Command: 'ls' })} setDoc={vi.fn()} />);

    editNode({ Command: 'ls', ProjectName: '' });

    expect(shownParams()).toEqual({ Command: 'ls', ProjectName: PROJECT });
  });

  it('fills a ProjectName seeded empty by a new node', () => {
    render(<DetailsViewDocument doc={makeWorkflowDoc({ Command: 'ls' })} setDoc={vi.fn()} />);

    act(() => graphProps.onAddNode());
    act(() => graphProps.onSetNode('node2', {
      ...graphProps.nodes.node2,
      Execution: { input_parameters: { ProjectName: '' } },
    }));

    expect(graphProps.nodes.node2.Execution.input_parameters.ProjectName).toBe(PROJECT);
  });

  it('fills an empty ProjectName that was already in the document', () => {
    render(<DetailsViewDocument doc={makeWorkflowDoc({ ProjectName: '' })} setDoc={vi.fn()} />);

    // The load path does not fill, so it is still empty on screen.
    expect(shownParams().ProjectName).toBe('');

    editNode({ ProjectName: '', Command: 'ls' });

    expect(shownParams().ProjectName).toBe(PROJECT);
  });

  it('matches the field name in any casing', () => {
    render(<DetailsViewDocument doc={makeWorkflowDoc({ Command: 'ls' })} setDoc={vi.fn()} />);

    editNode({ Command: 'ls', projectname: '' });

    expect(shownParams().projectname).toBe(PROJECT);
  });

  it('never blocks a user pointing the node at another project', () => {
    render(<DetailsViewDocument doc={makeWorkflowDoc({ ProjectName: '' })} setDoc={vi.fn()} />);

    editNode({ ProjectName: 'OTHER_PROJECT' });

    expect(shownParams().ProjectName).toBe('OTHER_PROJECT');
  });

  it('does not add a ProjectName field to a node that has none', () => {
    render(<DetailsViewDocument doc={makeWorkflowDoc({ Command: 'ls' })} setDoc={vi.fn()} />);

    editNode({ Command: 'pwd' });

    expect(shownParams()).toEqual({ Command: 'pwd' });
  });

  it('leaves the field empty when no project is selected', () => {
    useProjectStore.getState().selectProject(NO_PROJECT);
    render(<DetailsViewDocument doc={makeWorkflowDoc({ Command: 'ls' })} setDoc={vi.fn()} />);

    editNode({ Command: 'ls', ProjectName: '' });

    expect(shownParams().ProjectName).toBe('');
  });
});
