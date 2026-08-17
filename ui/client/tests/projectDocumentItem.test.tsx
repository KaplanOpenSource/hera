import { describe, it, expect, afterEach } from 'vitest';
import { render, screen, cleanup } from '@testing-library/react';
import { SimpleTreeView } from '@mui/x-tree-view';

const { ProjectDocumentItem } = await import('../src/components/project/ProjectDocumentItem');

afterEach(() => { cleanup(); });

const project = { name: 'TestProject', documents: [] } as any;
const document = {
  _id: { $oid: 'doc1' },
  desc: { datasourceName: 'MyDoc' },
  type: 'TestType',
  _cls: 'Metadata.Cache',
  resource: 'res',
  dataFormat: 'string',
  projectName: 'TestProject',
} as any;

const renderItem = () => render(
  <SimpleTreeView>
    <ProjectDocumentItem project={project} document={document} />
  </SimpleTreeView>
);

describe('ProjectDocumentItem', () => {
  it('renders the document name', () => {
    renderItem();
    expect(screen.getByText('MyDoc')).toBeTruthy();
  });

  it('shows a resource preview when the resource is a string', () => {
    renderItem();
    expect(screen.getByText(/resource: res/)).toBeTruthy();
  });

  it('does not render a delete button (delete lives in the document view)', () => {
    renderItem();
    expect(screen.queryByLabelText('Delete Document')).toBeNull();
  });
});
