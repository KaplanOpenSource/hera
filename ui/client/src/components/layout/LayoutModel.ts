import { Actions, DockLocation, IJsonModel, IJsonTabNode, Model, TabNode, TabSetNode } from 'flexlayout-react';
import { ProjectObj } from '../../objects/ProjectObj';
import { canvasTab } from './tabs/CanvasTab';
import { detailsTab } from './tabs/DetailsTab';
import { LayoutTab } from './tabs/LayoutTab';
import { outputTab } from './tabs/OutputTab';
import { previewTab } from './tabs/PreviewTab';
import { treeTab } from './tabs/TreeTab';
import { tabForId } from './tabs/allTabs';

// Tabset ids, and the id of the single tree tab.
export const TREE_TAB_ID = treeTab.tabId;
const TREE_TABSET_ID = 'tree-tabset';
const DETAILS_TABSET_ID = 'details-tabset';

const GLOBAL_LAYOUT_CONFIG = {
  tabEnableClose: true,
  tabEnableRename: true,
  tabEnableDrag: true,
  tabSetEnableMaximize: true,
  tabSetEnableClose: true,
  tabSetEnableDeleteWhenEmpty: true,
  rootOrientationVertical: false,
};

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
            ? [{ type: 'tabset', id: TREE_TABSET_ID, weight: 25, enableClose: false, enableDrag: false, enableDrop: false, enableTabStrip: false, children: [treeTab.node()] }]
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

  // All open tabs of one kind.
  private tabsOfKind(kind: LayoutTab): TabNode[] {
    const tabs: TabNode[] = [];
    this._model.visitNodes((node) => {
      if (node.getType() === 'tab' && kind.owns(node.getId())) {
        tabs.push(node as TabNode);
      }
    });
    return tabs;
  }

  // The document a tab shows, as its own kind reads it off the tab's config.
  private docIdOfNode(node: TabNode): string | undefined {
    return tabForId(node.getId())?.docIdOf(node.getConfig());
  }

  // The document a tab shows, or undefined if it shows no document.
  docIdOfTab(tabId: string): string | undefined {
    const tab = this.getTab(tabId);
    return tab ? this.docIdOfNode(tab) : undefined;
  }

  // Whether a tab other than the given one still shows the same document.
  hasOtherTabForDoc(docid: string, exceptTabId: string): boolean {
    let found = false;
    this._model.visitNodes((node) => {
      if (node.getType() !== 'tab' || node.getId() === exceptTabId) return;
      if (this.docIdOfNode(node as TabNode) === docid) {
        found = true;
      }
    });
    return found;
  }

  // Open the details tab for an item, or focus it if it is already open.
  openOrFocusDetailsTab(showItemId: string, project: ProjectObj): void {
    const detailsId = detailsTab.id(showItemId);
    if (this._model.getNodeById(detailsId)) {
      this._model.doAction(Actions.selectTab(detailsId));
    } else {
      this._model.doAction(Actions.addTab(detailsTab.node(showItemId, project), DETAILS_TABSET_ID, DockLocation.CENTER, -1));
    }
  }

  // Open a workflow's canvas below the details panel, or focus it if it is open.
  // Later canvases join the first one's tabset, so they don't each split the row.
  openOrFocusCanvasTab(docid: string, docName: string): void {
    const canvasId = canvasTab.id(docid);
    if (this._model.getNodeById(canvasId)) {
      this._model.doAction(Actions.selectTab(canvasId));
      return;
    }
    const tab = canvasTab.node(docid, docName);
    const openCanvas = this.tabsOfKind(canvasTab)[0];
    if (openCanvas) {
      this._model.doAction(Actions.addTab(tab, openCanvas.getParent()!.getId(), DockLocation.CENTER, -1));
    } else {
      this._model.doAction(Actions.addTab(tab, DETAILS_TABSET_ID, DockLocation.BOTTOM, -1));
    }
  }

  // Open a workflow's run output to the right of the canvas, or focus it if it is
  // open. Later outputs join the first one's tabset, like the canvas tabs do. With
  // no canvas open the output falls back to below the details panel.
  openOrFocusOutputTab(workflowName: string): void {
    const outputId = outputTab.id(workflowName);
    if (this._model.getNodeById(outputId)) {
      this._model.doAction(Actions.selectTab(outputId));
      return;
    }
    const tab = outputTab.node(workflowName);
    const openOutput = this.tabsOfKind(outputTab)[0];
    const openCanvas = this.tabsOfKind(canvasTab)[0];
    if (openOutput) {
      this._model.doAction(Actions.addTab(tab, openOutput.getParent()!.getId(), DockLocation.CENTER, -1));
    } else if (openCanvas) {
      this._model.doAction(Actions.addTab(tab, openCanvas.getParent()!.getId(), DockLocation.RIGHT, -1));
    } else {
      this._model.doAction(Actions.addTab(tab, DETAILS_TABSET_ID, DockLocation.BOTTOM, -1));
    }
  }

  // The details tab the user is looking at, or undefined when none is open. It is
  // the tab its tabset has selected, which also holds before the first layout pass.
  activeDetailsTab(): TabNode | undefined {
    let active: TabNode | undefined;
    this._model.visitNodes((node) => {
      if (node.getType() !== 'tabset') return;
      const selected = (node as TabSetNode).getSelectedNode();
      if (selected && detailsTab.owns(selected.getId())) {
        active = selected as TabNode;
      }
    });
    return active;
  }

  // Close a document's canvas tab, if it has one open.
  closeCanvasTab(docid: string): void {
    const tab = this.getTab(canvasTab.id(docid));
    if (tab) {
      this._model.doAction(Actions.deleteTab(tab.getId()));
    }
  }

  // Show or hide the tree panel.
  setTreeVisible(visible: boolean): void {
    const treeNode = this._model.getNodeById(TREE_TAB_ID);
    if (!visible && treeNode) {
      this._model.doAction(Actions.deleteTab(TREE_TAB_ID));
    } else if (visible && !treeNode) {
      this._model.doAction(Actions.addTab(treeTab.node(), DETAILS_TABSET_ID, DockLocation.LEFT, -1));
    }
  }

  // Replace the preview pane: clear any existing preview tab, then open one for
  // the given document if both a docid and name are supplied.
  setPreview(docid?: string, docName?: string): void {
    for (const t of this.tabsOfKind(previewTab)) {
      this._model.doAction(Actions.deleteTab(t.getId()));
    }
    if (docid && docName !== undefined) {
      this._model.doAction(Actions.addTab(previewTab.node(docid, docName), DETAILS_TABSET_ID, DockLocation.BOTTOM, -1));
    }
  }

  // A fresh model in the reset arrangement that keeps the open details tabs,
  // preserving which one is active.
  resetKeepingDetails(treeVisible: boolean, activeShowItemId: string | undefined): LayoutModel {
    const detailsTabs = this.tabsOfKind(detailsTab).map(t => t.toJson());
    const selectedIndex = detailsTabs.findIndex(t => t.id === detailsTab.id(activeShowItemId ?? ''));
    return LayoutModel.create(treeVisible, detailsTabs, selectedIndex);
  }

  // Close details/preview tabs whose document no longer exists (e.g. after it was deleted).
  closeMissingDocuments(project: ProjectObj): void {
    for (const kind of [detailsTab, previewTab, canvasTab]) {
      for (const t of this.tabsOfKind(kind)) {
        const oid = kind.docIdOf(t.getConfig());
        if (oid && !project.documentIds.has(oid)) {
          this._model.doAction(Actions.deleteTab(t.getId()));
        }
      }
    }
  }

  // Rename open details tabs whose computed name has drifted from the project. A
  // tab opened for a just-created document (e.g. a new notebook) is named before
  // that document has loaded, so detailsTabName falls back to the project-config
  // name; once the document arrives, rename the tab to its real name.
  syncTabNames(project: ProjectObj): void {
    for (const t of this.tabsOfKind(detailsTab)) {
      const showItemId = t.getConfig()?.showItemId as string | undefined;
      if (!showItemId) continue;
      const name = detailsTab.tabName(showItemId, project);
      if (name !== t.getName()) {
        this._model.doAction(Actions.renameTab(t.getId(), name));
      }
    }
  }
}
