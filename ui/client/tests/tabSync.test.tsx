import { describe, it, expect } from 'vitest';
import { IJsonTabNode } from 'flexlayout-react';
import { LayoutModel } from '../src/components/layout/LayoutModel';
import { ProjectObj } from '../src/objects/ProjectObj';
import { idDocId } from '../src/shared/idDocId';
import { ProjectDocument } from '../src/shared/types';

const notebookDoc = (oid: string, notebookName: string): ProjectDocument => ({
  _cls: 'Metadata.Cache',
  _id: { $oid: oid },
  projectName: 'TestProject',
  desc: { datasourceName: notebookName },
  type: 'notebook',
  resource: `/notebooks/${notebookName}.ipynb`,
  dataFormat: 'JSON_dict',
});

const projectWith = (docs: ProjectDocument[]): ProjectObj => {
  return new ProjectObj({ name: 'TestProject', documents: docs });
};

const detailsTab = (oid: string, name: string): IJsonTabNode => {
  const showItemId = idDocId(oid);
  return { type: 'tab', id: `details:${showItemId}`, name, component: 'details', config: { showItemId } };
};

describe('LayoutModel.syncTabNames', () => {
  it('renames a new notebook tab once its document loads', () => {
    const oid = 'nb1';
    const tabId = `details:${idDocId(oid)}`;
    // The tab was opened before the notebook document loaded, so its name fell
    // back to the project-config name.
    const staleName = 'TestProject config';
    const layout = LayoutModel.create(true, [detailsTab(oid, staleName)]);

    // While the document is still missing from the project, the fallback name
    // matches, so nothing is renamed.
    layout.syncTabNames(projectWith([]));
    expect(layout.getTab(tabId)?.getName()).toBe(staleName);

    // Once the notebook document arrives, the tab takes its real name.
    layout.syncTabNames(projectWith([notebookDoc(oid, 'myNotebook')]));
    expect(layout.getTab(tabId)?.getName()).toBe('myNotebook');
  });

  it('leaves an already-correct tab name unchanged', () => {
    const oid = 'nb2';
    const tabId = `details:${idDocId(oid)}`;
    const layout = LayoutModel.create(true, [detailsTab(oid, 'analysis')]);
    layout.syncTabNames(projectWith([notebookDoc(oid, 'analysis')]));
    expect(layout.getTab(tabId)?.getName()).toBe('analysis');
  });
});

describe('LayoutModel.closeMissingDocuments', () => {
  it('closes a details tab whose document no longer exists', () => {
    const oid = 'gone';
    const tabId = `details:${idDocId(oid)}`;
    const layout = LayoutModel.create(true, [detailsTab(oid, 'Old Doc')]);
    layout.closeMissingDocuments(projectWith([]));
    expect(layout.getTab(tabId)).toBeUndefined();
  });

  it('keeps a details tab whose document still exists', () => {
    const oid = 'keep';
    const tabId = `details:${idDocId(oid)}`;
    const layout = LayoutModel.create(true, [detailsTab(oid, 'keeper')]);
    layout.closeMissingDocuments(projectWith([notebookDoc(oid, 'keeper')]));
    expect(layout.getTab(tabId)).toBeDefined();
  });

  it('closes a preview tab whose document no longer exists', () => {
    const oid = 'previewgone';
    const previewTab: IJsonTabNode = {
      type: 'tab',
      id: `preview:${oid}`,
      name: 'Preview: Old',
      component: 'preview',
      config: { docid: oid },
    };
    const layout = LayoutModel.create(true, [previewTab]);
    layout.closeMissingDocuments(projectWith([]));
    expect(layout.getTab(`preview:${oid}`)).toBeUndefined();
  });
});
