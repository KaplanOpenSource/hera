import { Actions, DockLocation, IJsonModel, IJsonTabNode, Model, TabNode } from 'flexlayout-react';
import { ProjectObj } from '../../objects/ProjectObj';
import { idFromDocId } from '../../shared/idDocId';
import { tabKindClassName } from '../../shared/tabKind';
import { detailsTabName } from '../details/DetailsViewPanel';
import { LayoutComponent } from './LayoutPanel';

// Tabset ids, the single tree tab id, and the id prefixes for the per-item tabs.
export const TREE_TAB_ID = 'tree';
const TREE_TABSET_ID = 'tree-tabset';
const DETAILS_TABSET_ID = 'details-tabset';
export const DETAILS_TAB_PREFIX = 'details:';
const PREVIEW_TAB_PREFIX = 'preview:';

const GLOBAL_LAYOUT_CONFIG = {
  tabEnableClose: true,
  tabEnableRename: true,
  tabEnableDrag: true,
  tabSetEnableMaximize: true,
  tabSetEnableClose: true,
  tabSetEnableDeleteWhenEmpty: true,
  rootOrientationVertical: false,
};

const TREE_TAB: IJsonTabNode = { type: 'tab', id: TREE_TAB_ID, name: 'Project', component: LayoutComponent.Tree };

// Tab node for a tree item's details view.
const makeDetailsTab = (showItemId: string, project: ProjectObj): IJsonTabNode => ({
  type: 'tab',
  id: `${DETAILS_TAB_PREFIX}${showItemId}`,
  name: detailsTabName(showItemId, project),
  className: tabKindClassName(showItemId, project),
  component: LayoutComponent.Details,
  config: { showItemId },
});

// Tab node for a document's preview pane.
const makePreviewTab = (docid: string, docName: string): IJsonTabNode => ({
  type: 'tab',
  id: `${PREVIEW_TAB_PREFIX}${docid}`,
  name: `Preview: ${docName}`,
  component: LayoutComponent.Preview,
  config: { docid },
});

// Wraps a flexlayout Model, exposing the layout operations this app performs on
// it. The raw model is reached via `.model` to hand to <Layout model={...}>.
export class LayoutModel {
  private constructor(private readonly _model: Model) {}

  // Builds a fresh model in the original arrangement: tree panel (25%) on the
  // left, details panel (75%) on the right. Any supplied details tabs are placed
  // in the details panel so a layout reset doesn't lose the user's open documents.
  static create(treeVisible: boolean, detailsTabs: IJsonTabNode[] = [], selectedIndex = -1): LayoutModel {
    const detailsTabset = {
      type: 'tabset',
      id: DETAILS_TABSET_ID,
      weight: 75,
      enableDeleteWhenEmpty: false,
      children: detailsTabs,
      ...(selectedIndex >= 0 ? { selected: selectedIndex } : {}),
    };
    const layout: IJsonModel = {
      global: GLOBAL_LAYOUT_CONFIG,
      layout: {
        type: 'row',
        children: [
          ...(treeVisible
            ? [{ type: 'tabset', id: TREE_TABSET_ID, weight: 25, children: [TREE_TAB] }]
            : []),
          detailsTabset,
        ],
      },
    };
    return new LayoutModel(Model.fromJson(layout));
  }

  // The underlying flexlayout model, to pass to <Layout model={...}>.
  get model(): Model {
    return this._model;
  }

  // The tab node with the given id, or undefined if there is no such tab.
  getTab(id: string): TabNode | undefined {
    const node = this._model.getNodeById(id);
    return node?.getType() === 'tab' ? (node as TabNode) : undefined;
  }

  // All open tabs whose id starts with the given prefix (details or preview tabs).
  private tabsWithPrefix(prefix: string): TabNode[] {
    const tabs: TabNode[] = [];
    this._model.visitNodes((node) => {
      if (node.getType() === 'tab' && node.getId().startsWith(prefix)) {
        tabs.push(node as TabNode);
      }
    });
    return tabs;
  }

  // Open the details tab for an item, or focus it if it is already open.
  openOrFocusDetailsTab(showItemId: string, project: ProjectObj): void {
    const detailsId = `${DETAILS_TAB_PREFIX}${showItemId}`;
    if (this._model.getNodeById(detailsId)) {
      this._model.doAction(Actions.selectTab(detailsId));
    } else {
      this._model.doAction(Actions.addTab(makeDetailsTab(showItemId, project), DETAILS_TABSET_ID, DockLocation.CENTER, -1));
    }
  }

  // Show or hide the tree panel.
  setTreeVisible(visible: boolean): void {
    const treeNode = this._model.getNodeById(TREE_TAB_ID);
    if (!visible && treeNode) {
      this._model.doAction(Actions.deleteTab(TREE_TAB_ID));
    } else if (visible && !treeNode) {
      this._model.doAction(Actions.addTab(TREE_TAB, DETAILS_TABSET_ID, DockLocation.LEFT, -1));
    }
  }

  // Replace the preview pane: clear any existing preview tab, then open one for
  // the given document if both a docid and name are supplied.
  setPreview(docid?: string, docName?: string): void {
    for (const t of this.tabsWithPrefix(PREVIEW_TAB_PREFIX)) {
      this._model.doAction(Actions.deleteTab(t.getId()));
    }
    if (docid && docName !== undefined) {
      this._model.doAction(Actions.addTab(makePreviewTab(docid, docName), DETAILS_TABSET_ID, DockLocation.BOTTOM, -1));
    }
  }

  // A fresh model in the reset arrangement that keeps the open details tabs,
  // preserving which one is active.
  resetKeepingDetails(treeVisible: boolean, activeShowItemId: string | undefined): LayoutModel {
    const detailsTabs = this.tabsWithPrefix(DETAILS_TAB_PREFIX).map(t => t.toJson());
    const selectedIndex = detailsTabs.findIndex(t => t.id === `${DETAILS_TAB_PREFIX}${activeShowItemId}`);
    return LayoutModel.create(treeVisible, detailsTabs, selectedIndex);
  }

  // Close details/preview tabs whose document no longer exists (e.g. after it was deleted).
  closeMissingDocuments(project: ProjectObj): void {
    for (const t of this.tabsWithPrefix(DETAILS_TAB_PREFIX)) {
      const oid = idFromDocId(t.getConfig()?.showItemId ?? '');
      if (oid && !project.documentIds.has(oid)) {
        this._model.doAction(Actions.deleteTab(t.getId()));
      }
    }
    for (const t of this.tabsWithPrefix(PREVIEW_TAB_PREFIX)) {
      const oid = t.getConfig()?.docid as string | undefined;
      if (oid && !project.documentIds.has(oid)) {
        this._model.doAction(Actions.deleteTab(t.getId()));
      }
    }
  }

  // Rename open details tabs whose computed name has drifted from the project. A
  // tab opened for a just-created document (e.g. a new notebook) is named before
  // that document has loaded, so detailsTabName falls back to the project-config
  // name; once the document arrives, rename the tab to its real name.
  syncTabNames(project: ProjectObj): void {
    for (const t of this.tabsWithPrefix(DETAILS_TAB_PREFIX)) {
      const showItemId = t.getConfig()?.showItemId as string | undefined;
      if (!showItemId) continue;
      const name = detailsTabName(showItemId, project);
      if (name !== t.getName()) {
        this._model.doAction(Actions.renameTab(t.getId(), name));
      }
    }
  }
}
