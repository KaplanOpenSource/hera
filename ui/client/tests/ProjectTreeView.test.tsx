import { describe, it, expect, vi, afterEach } from 'vitest';
import { render, screen, cleanup } from '@testing-library/react';
import { ProjectTreeView } from '../src/components/project/ProjectTreeView';
import { ProjectObj } from '../src/objects/ProjectObj';
import { ProjectEntire } from '../src/shared/types';

vi.mock('../src/io/FetchProjects', () => ({
  fetchProjectDetails: vi.fn(),
}));
vi.mock('../src/io/snackbar', () => ({
  pushRunning: () => 'mock-key',
  pushError: () => 'mock-key',
  dismiss: () => {},
}));
vi.mock('../src/components/project/AddDocumentButton', () => ({
  AddDocumentButton: () => null,
}));
vi.mock('../src/components/project/ProjectViewSettingsButton', () => ({
  ProjectViewSettingsButton: () => null,
}));
vi.mock('../src/components/project/RepoTreeWhole', () => ({
  RepoTreeWhole: () => null,
}));

const makeProject = (documents: ProjectEntire['documents'] = []): ProjectObj => {
  return new ProjectObj({
    name: 'TestProject',
    documents: [
      {
        _cls: 'Config',
        _id: { $oid: 'config-id' },
        projectName: 'TestProject',
        desc: { filesDirectory: '/tmp/test' },
        type: 'TestProject__config__',
        resource: '',
        dataFormat: 'json',
      },
      ...documents,
    ],
  });
};

const renderTree = (project: ProjectObj) => {
  return render(
    <ProjectTreeView
      project={project}
      setSelectedItemIds={vi.fn()}
    />
  );
};

describe('ProjectTreeView', () => {
  afterEach(() => {
    cleanup();
  });

  it('renders the project name', () => {
    renderTree(makeProject());
    expect(screen.getByText('Project TestProject')).toBeDefined();
  });

  it('renders the files directory', () => {
    renderTree(makeProject());
    expect(screen.getByText('/tmp/test')).toBeDefined();
  });

  it('renders document names when documents exist', () => {
    const project = makeProject([
      {
        _cls: 'Document',
        _id: { $oid: 'doc-1' },
        projectName: 'TestProject',
        desc: { datasourceName: 'MyDataSource', toolkit: 'toolkitA' },
        type: 'typeA',
        resource: 'some resource',
        dataFormat: 'json',
      },
    ]);
    renderTree(project);
    expect(screen.getByText('MyDataSource')).toBeDefined();
  });

  it('groups documents by toolkit when multiple toolkits exist', () => {
    const project = makeProject([
      {
        _cls: 'Document',
        _id: { $oid: 'doc-1' },
        projectName: 'TestProject',
        desc: { datasourceName: 'Doc1', toolkit: 'toolkitA' },
        type: 'typeA',
        resource: '',
        dataFormat: 'json',
      },
      {
        _cls: 'Document',
        _id: { $oid: 'doc-2' },
        projectName: 'TestProject',
        desc: { datasourceName: 'Doc2', toolkit: 'toolkitA' },
        type: 'typeA',
        resource: '',
        dataFormat: 'json',
      },
      {
        _cls: 'Document',
        _id: { $oid: 'doc-3' },
        projectName: 'TestProject',
        desc: { datasourceName: 'Doc3', toolkit: 'toolkitB' },
        type: 'typeA',
        resource: '',
        dataFormat: 'json',
      },
      {
        _cls: 'Document',
        _id: { $oid: 'doc-4' },
        projectName: 'TestProject',
        desc: { datasourceName: 'Doc4', toolkit: 'toolkitB' },
        type: 'typeA',
        resource: '',
        dataFormat: 'json',
      },
    ]);
    renderTree(project);
    expect(screen.getByText('toolkitA')).toBeDefined();
    expect(screen.getByText('toolkitB')).toBeDefined();
  });

  it('renders documents flat when no differentiating desc fields exist', () => {
    const project = makeProject([
      {
        _cls: 'Document',
        _id: { $oid: 'doc-1' },
        projectName: 'TestProject',
        desc: { datasourceName: 'Doc1' },
        type: 'typeA',
        resource: '',
        dataFormat: 'json',
      },
      {
        _cls: 'Document',
        _id: { $oid: 'doc-2' },
        projectName: 'TestProject',
        desc: { datasourceName: 'Doc2' },
        type: 'typeA',
        resource: '',
        dataFormat: 'json',
      },
    ]);
    renderTree(project);
    expect(screen.getByText('Doc1')).toBeDefined();
    expect(screen.getByText('Doc2')).toBeDefined();
  });
});
