import { describe, it, expect } from 'vitest';
import { IJsonTabNode, Model, TabNode } from 'flexlayout-react';
import {
  createLayoutModel,
  syncDetailsTabNames,
  closeTabsForMissingDocuments,
} from '../src/components/layout/ProjectLayout';
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

const tabName = (model: Model, id: string): string | undefined => {
  return (model.getNodeById(id) as TabNode | undefined)?.getName();
};

describe('syncDetailsTabNames', () => {
  it('renames a new notebook tab once its document loads', () => {
    const oid = 'nb1';
    const tabId = `details:${idDocId(oid)}`;
    // The tab was opened before the notebook document loaded, so its name fell
    // back to the project-config name.
    const staleName = 'TestProject config';
    const model = createLayoutModel(true, [detailsTab(oid, staleName)]);

    // While the document is still missing from the project, the fallback name
    // matches, so nothing is renamed.
    syncDetailsTabNames(model, projectWith([]));
    expect(tabName(model, tabId)).toBe(staleName);

    // Once the notebook document arrives, the tab takes its real name.
    syncDetailsTabNames(model, projectWith([notebookDoc(oid, 'myNotebook')]));
    expect(tabName(model, tabId)).toBe('myNotebook');
  });

  it('leaves an already-correct tab name unchanged', () => {
    const oid = 'nb2';
    const tabId = `details:${idDocId(oid)}`;
    const model = createLayoutModel(true, [detailsTab(oid, 'analysis')]);
    syncDetailsTabNames(model, projectWith([notebookDoc(oid, 'analysis')]));
    expect(tabName(model, tabId)).toBe('analysis');
  });
});

describe('closeTabsForMissingDocuments', () => {
  it('closes a details tab whose document no longer exists', () => {
    const oid = 'gone';
    const tabId = `details:${idDocId(oid)}`;
    const model = createLayoutModel(true, [detailsTab(oid, 'Old Doc')]);
    closeTabsForMissingDocuments(model, projectWith([]));
    expect(model.getNodeById(tabId)).toBeUndefined();
  });

  it('keeps a details tab whose document still exists', () => {
    const oid = 'keep';
    const tabId = `details:${idDocId(oid)}`;
    const model = createLayoutModel(true, [detailsTab(oid, 'keeper')]);
    closeTabsForMissingDocuments(model, projectWith([notebookDoc(oid, 'keeper')]));
    expect(model.getNodeById(tabId)).toBeDefined();
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
    const model = createLayoutModel(true, [previewTab]);
    closeTabsForMissingDocuments(model, projectWith([]));
    expect(model.getNodeById(`preview:${oid}`)).toBeUndefined();
  });
});
