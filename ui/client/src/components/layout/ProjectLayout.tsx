import { Box, Paper } from '@mui/material';
import { Action, Actions, DockLocation, IJsonModel, IJsonTabNode, ITabRenderValues, Layout, Model, TabNode } from 'flexlayout-react';
import 'flexlayout-react/style/light.css';
import { useCallback, useEffect, useState } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
import { ProjectObj } from '../../objects/ProjectObj';
import { CENTRAL_REPO_FOLDER_ID, idFromDocId, idFromRepoId, isSplitId, normalizeSplitId } from '../../shared/idDocId';
import { classifyTab, tabKindClassName } from '../../shared/tabKind';
import { TAB_KIND_STYLES, tabKindCss } from '../../shared/tabKindConfig';
import { DetailsViewPanel, detailsTabName } from '../details/DetailsViewPanel';
import { hasPreview, PreviewPanel } from '../details/PreviewPanel';
import { ProjectTreeView } from '../project/ProjectTreeView';

const GLOBAL_LAYOUT_CONFIG = {
  tabEnableClose: true,
  tabEnableRename: true,
  tabEnableDrag: true,
  tabSetEnableMaximize: true,
  tabSetEnableClose: true,
  tabSetEnableDeleteWhenEmpty: true,
  rootOrientationVertical: false,
};

const TREE_TAB: IJsonTabNode = { type: 'tab', id: 'tree', name: 'Project', component: 'tree' };

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
    id: 'details-tabset',
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
          ? [{ type: 'tabset', id: 'tree-tabset', weight: 25, children: [TREE_TAB] }]
          : []),
        detailsTabset,
      ],
    },
  };
  return Model.fromJson(layout);
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
    docId ? [`document_${docId}`] : []
  );
  const [activeShowItemId, setActiveShowItemId] = useState<string | undefined>(
    docId ? `document_${docId}` : undefined
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
    const detailsTabs: IJsonTabNode[] = [];
    model.visitNodes((node) => {
      if (node.getType() === 'tab' && node.getId().startsWith('details:')) {
        detailsTabs.push((node as TabNode).toJson());
      }
    });
    const selectedIndex = detailsTabs.findIndex(t => t.id === `details:${activeShowItemId}`);
    setModel(createLayoutModel(!treeCollapsed, detailsTabs, selectedIndex));
  }, [resetSignal]);

  useEffect(() => {
    const treeNode = model.getNodeById('tree');
    if (treeCollapsed && treeNode) {
      model.doAction(Actions.deleteTab('tree'));
    } else if (!treeCollapsed && !treeNode) {
      model.doAction(Actions.addTab(
        { type: 'tab', id: 'tree', name: 'Project', component: 'tree' },
        'details-tabset',
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
      : 'config';

    const detailsId = `details:${showItemId}`;
    const existing = model.getNodeById(detailsId);
    if (existing) {
      model.doAction(Actions.selectTab(detailsId));
    } else {
      model.doAction(Actions.addTab(
        {
          type: 'tab',
          id: detailsId,
          name: detailsTabName(showItemId, project),
          className: tabKindClassName(showItemId, project),
          component: 'details',
          config: { showItemId },
        },
        'details-tabset',
        DockLocation.CENTER,
        -1,
      ));
    }
    setActiveShowItemId(showItemId);
  }, [selectedItemsIds[0], model]);

  useEffect(() => {
    const existingPreview: string[] = [];
    model.visitNodes((node) => {
      if (node.getType() === 'tab' && node.getId().startsWith('preview:')) {
        existingPreview.push(node.getId());
      }
    });
    existingPreview.forEach(id => model.doAction(Actions.deleteTab(id)));

    if (previewAvailable && activeDocId) {
      model.doAction(Actions.addTab(
        {
          type: 'tab',
          id: `preview:${activeDocId}`,
          name: `Preview: ${activeDoc!.name}`,
          component: 'preview',
          config: { docid: activeDocId },
        },
        'details-tabset',
        DockLocation.BOTTOM,
        -1,
      ));
    }
  }, [activeShowItemId, previewAvailable, activeDocId, model]);

  useEffect(() => {
    if (docId && project?.documentIds.has(docId)) {
      setSelectedItemIds([`document_${docId}`]);
    } else {
      setSelectedItemIds([]);
    }
  }, [project?.name]);

  const handleSetSelectedItemIds = useCallback((ids: string[]) => {
    setSelectedItemIds(ids);
    const selectedId = ids[0];
    const oid = selectedId?.startsWith('document_') ? selectedId.slice('document_'.length) : undefined;
    const basePath = '/' + encodeURIComponent(project?.name ?? '');
    const newPath = oid ? `${basePath}/${oid}` : basePath;
    if (location.pathname !== newPath) {
      navigate(newPath, { replace: true });
    }
  }, [project?.name, navigate]);

  const handleAction = useCallback((action: Action) => {
    if (action.type === Actions.SELECT_TAB) {
      const tabId = action.data.tabNode as string;
      if (tabId?.startsWith('details:')) {
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
      case 'tree':
        return (
          <Paper sx={{ p: 2, height: '100%', overflow: 'auto' }}>
            <ProjectTreeView
              project={project}
              selectedItemsIds={selectedItemsIds}
              setSelectedItemIds={handleSetSelectedItemIds}
            />
          </Paper>
        );
      case 'details':
        return (
          <Paper sx={{ height: '100%', overflow: 'hidden' }}>
            <DetailsViewPanel
              project={project}
              showItemId={config?.showItemId}
            />
          </Paper>
        );
      case 'preview':
        return <PreviewPanel docid={config?.docid} />;
      default:
        return null;
    }
  };

  return (
    <Box sx={{ position: 'relative', flex: 1, height: '100%' }}>
      <style>{tabKindCss}</style>
      <Layout
        model={model}
        factory={factory}
        onRenderTab={onRenderTab}
        onAction={handleAction}
      />
    </Box>
  );
};
