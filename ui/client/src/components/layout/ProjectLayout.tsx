import { Paper } from '@mui/material';
import { Action, Actions, DockLocation, IJsonModel, IJsonTabNode, ITabRenderValues, Layout, Model, TabNode } from 'flexlayout-react';
import 'flexlayout-react/style/light.css';
import { useCallback, useEffect, useState } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
import { ProjectObj } from '../../objects/ProjectObj';
import { CENTRAL_REPO_FOLDER_ID, idDocId, idFromDocId, idFromRepoId, isSplitId, normalizeSplitId } from '../../shared/idDocId';
import { classifyTab, tabKindClassName } from '../../shared/tabKind';
import { TAB_KIND_STYLES } from '../../shared/tabKindConfig';
import { DetailsViewPanel, detailsTabName } from '../details/DetailsViewPanel';
import { hasPreview, PreviewPanel } from '../details/PreviewPanel';
import { ProjectTreeView } from '../project/ProjectTreeView';

// The component identifiers flexlayout stores on each tab and passes to the factory.
enum LayoutComponent {
  Tree = 'tree',
  Details = 'details',
  Preview = 'preview',
}

// Tabset ids, the single tree tab id, and the id prefixes for the per-item tabs.
const TREE_TAB_ID = 'tree';
const TREE_TABSET_ID = 'tree-tabset';
const DETAILS_TABSET_ID = 'details-tabset';
const DETAILS_TAB_PREFIX = 'details:';
const PREVIEW_TAB_PREFIX = 'preview:';
// Fallback "item" shown when the tree selection isn't a document/repo/split.
const CONFIG_ITEM_ID = 'config';

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

// Builds a fresh layout model in the original arrangement: tree panel (25%) on the
// left, details panel (75%) on the right. Any open details tabs are re-placed into
// the details panel so a layout reset doesn't lose the user's open documents.
export const createLayoutModel = (
  treeVisible: boolean,
  detailsTabs: IJsonTabNode[] = [],
  selectedIndex = -1,
): Model => {
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
  return Model.fromJson(layout);
};

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

// All open tabs whose id starts with the given prefix (e.g. all details or preview tabs).
const tabsWithPrefix = (model: Model, prefix: string): TabNode[] => {
  const tabs: TabNode[] = [];
  model.visitNodes((node) => {
    if (node.getType() === 'tab' && node.getId().startsWith(prefix)) {
      tabs.push(node as TabNode);
    }
  });
  return tabs;
};

export const ProjectLayout = ({
  project,
  treeCollapsed,
  resetSignal,
}: {
  project: ProjectObj,
  treeCollapsed: boolean,
  resetSignal: number,
}) => {
  const { docId } = useParams<{ docId: string }>();
  const navigate = useNavigate();
  const [selectedItemsIds, setSelectedItemIds] = useState<string[]>(
    docId ? [idDocId(docId)] : []
  );
  const [activeShowItemId, setActiveShowItemId] = useState<string | undefined>(
    docId ? idDocId(docId) : undefined
  );

  const [model, setModel] = useState(() => createLayoutModel(!treeCollapsed));

  const activeDocId = activeShowItemId ? idFromDocId(activeShowItemId) : undefined;
  const activeDoc = activeDocId
    ? project.allDocuments.find(d => d.docid === activeDocId)
    : undefined;
  const previewAvailable = activeDoc ? hasPreview(activeDoc.data) : false;

  // When the reset button in the header is pressed, rebuild the layout in its
  // original left/right arrangement while keeping the currently open details tabs.
  useEffect(() => {
    if (resetSignal === 0) return;
    const detailsTabs = tabsWithPrefix(model, DETAILS_TAB_PREFIX).map(t => t.toJson());
    const selectedIndex = detailsTabs.findIndex(t => t.id === `${DETAILS_TAB_PREFIX}${activeShowItemId}`);
    setModel(createLayoutModel(!treeCollapsed, detailsTabs, selectedIndex));
  }, [resetSignal]);

  useEffect(() => {
    const treeNode = model.getNodeById(TREE_TAB_ID);
    if (treeCollapsed && treeNode) {
      model.doAction(Actions.deleteTab(TREE_TAB_ID));
    } else if (!treeCollapsed && !treeNode) {
      model.doAction(Actions.addTab(
        TREE_TAB,
        DETAILS_TABSET_ID,
        DockLocation.LEFT,
        -1,
      ));
    }
  }, [treeCollapsed, model]);

  useEffect(() => {
    const rawShowItemId = selectedItemsIds[0];
    if (!rawShowItemId) return;

    const isSpecific = rawShowItemId === CENTRAL_REPO_FOLDER_ID
      || !!idFromDocId(rawShowItemId)
      || !!idFromRepoId(rawShowItemId)
      || isSplitId(rawShowItemId);
    const showItemId = isSpecific
      ? (isSplitId(rawShowItemId) ? normalizeSplitId(rawShowItemId, project.documents) : rawShowItemId)
      : CONFIG_ITEM_ID;

    const detailsId = `${DETAILS_TAB_PREFIX}${showItemId}`;
    const existing = model.getNodeById(detailsId);
    if (existing) {
      model.doAction(Actions.selectTab(detailsId));
    } else {
      model.doAction(Actions.addTab(
        makeDetailsTab(showItemId, project),
        DETAILS_TABSET_ID,
        DockLocation.CENTER,
        -1,
      ));
    }
    setActiveShowItemId(showItemId);
  }, [selectedItemsIds[0], model]);

  useEffect(() => {
    tabsWithPrefix(model, PREVIEW_TAB_PREFIX)
      .forEach(t => model.doAction(Actions.deleteTab(t.getId())));

    if (previewAvailable && activeDocId) {
      model.doAction(Actions.addTab(
        makePreviewTab(activeDocId, activeDoc!.name),
        DETAILS_TABSET_ID,
        DockLocation.BOTTOM,
        -1,
      ));
    }
  }, [activeShowItemId, previewAvailable, activeDocId, model]);

  useEffect(() => {
    if (docId && project?.documentIds.has(docId)) {
      setSelectedItemIds([idDocId(docId)]);
    } else {
      setSelectedItemIds([]);
    }
  }, [project?.name]);

  const handleSetSelectedItemIds = useCallback((ids: string[]) => {
    setSelectedItemIds(ids);
    const selectedId = ids[0];
    const oid = selectedId ? idFromDocId(selectedId) : undefined;
    const basePath = '/' + encodeURIComponent(project?.name ?? '');
    const newPath = oid ? `${basePath}/${oid}` : basePath;
    if (location.pathname !== newPath) {
      navigate(newPath, { replace: true });
    }
  }, [project?.name, navigate]);

  const handleAction = useCallback((action: Action) => {
    if (action.type === Actions.SELECT_TAB) {
      const tabId = action.data.tabNode as string;
      if (tabId?.startsWith(DETAILS_TAB_PREFIX)) {
        const node = model.getNodeById(tabId) as TabNode | undefined;
        if (node) {
          setActiveShowItemId(node.getConfig()?.showItemId);
        }
      }
    }
    return action;
  }, [model]);

  const onRenderTab = useCallback((node: TabNode, renderValues: ITabRenderValues) => {
    const showItemId = node.getConfig()?.showItemId as string | undefined;
    if (!showItemId) return;
    const kind = classifyTab(showItemId, project);
    if (!kind) return;
    const { icon: Icon, color } = TAB_KIND_STYLES[kind];
    renderValues.leading = <Icon sx={{ fontSize: 16, color }} />;
  }, [project.allDocuments]);

  const factory = (node: TabNode) => {
    const component = node.getComponent();
    const config = node.getConfig();
    switch (component) {
      case LayoutComponent.Tree:
        return (
          <Paper sx={{ p: 2, height: '100%', overflow: 'auto' }}>
            <ProjectTreeView
              project={project}
              selectedItemsIds={selectedItemsIds}
              setSelectedItemIds={handleSetSelectedItemIds}
            />
          </Paper>
        );
      case LayoutComponent.Details:
        return (
          <Paper sx={{ height: '100%', overflow: 'hidden' }}>
            <DetailsViewPanel
              project={project}
              showItemId={config?.showItemId}
            />
          </Paper>
        );
      case LayoutComponent.Preview:
        return <PreviewPanel docid={config?.docid} />;
      default:
        return null;
    }
  };

  return (
    <Layout
      model={model}
      factory={factory}
      onRenderTab={onRenderTab}
      onAction={handleAction}
    />
  );
};
