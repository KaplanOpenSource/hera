import { describe, it, expect } from 'vitest';
import { DocumentObj } from '../src/objects/DocumentObj';
import { ProjectObj } from '../src/objects/ProjectObj';
import { classifyDocument, TabKind } from '../src/shared/tabKind';

const PROJECT = 'P';
const project = new ProjectObj({ name: PROJECT, documents: [] } as any);

const makeDoc = (data: any) => new DocumentObj(
  {
    _id: { $oid: 'x' },
    projectName: PROJECT,
    desc: {},
    _cls: 'Some.Cls',
    type: 't',
    resource: 'r',
    dataFormat: 'string',
    ...data,
  },
  project,
);

describe('classifyDocument', () => {
  it('classifies the project config document', () => {
    expect(classifyDocument(makeDoc({ type: `${PROJECT}__config__` }))).toBe(TabKind.ProjectConfig);
  });

  it('classifies a notebook document', () => {
    expect(classifyDocument(makeDoc({ type: 'notebook' }))).toBe(TabKind.Notebook);
  });

  it('classifies a workflow document by its type', () => {
    expect(classifyDocument(makeDoc({ type: 'hermesWorkflow' }))).toBe(TabKind.Workflow);
  });

  it('classifies an agent document (resource holds effects)', () => {
    expect(classifyDocument(makeDoc({ resource: { effects: {} } }))).toBe(TabKind.Agent);
  });

  it('classifies a plain document', () => {
    expect(classifyDocument(makeDoc({ type: 'ToolkitDataSource', resource: 'abc' }))).toBe(TabKind.Document);
  });

  it('prefers config over agent when the type is the config type', () => {
    const doc = makeDoc({ type: `${PROJECT}__config__`, resource: { effects: {} } });
    expect(classifyDocument(doc)).toBe(TabKind.ProjectConfig);
  });
});
