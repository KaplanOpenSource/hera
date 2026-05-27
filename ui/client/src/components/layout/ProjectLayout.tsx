import { Code, Description, FolderOpen, Science, Settings, Source } from '@mui/icons-material';
import { Box, Paper } from '@mui/material';
import { Action, Actions, DockLocation, ITabRenderValues, Layout, Model, TabNode } from 'flexlayout-react';
import 'flexlayout-react/style/light.css';
import { useCallback, useEffect, useState } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
import { ProjectObj } from '../../objects/ProjectObj';
import { CENTRAL_REPO_FOLDER_ID, idFromDocId, idFromRepoId } from '../../shared/idDocId';
import { classifyTab, TabKind } from '../../shared/tabKind';
import { DetailsViewPanel, detailsTabName } from '../details/DetailsViewPanel';
import { hasPreview, PreviewPanel } from '../details/PreviewPanel';
import { ProjectTreeView } from '../project/ProjectTreeView';

export const ProjectLayout = ({
  project,
  treeCollapsed,
}: {
  project: ProjectObj,
  treeCollapsed: boolean,
}) => {
  const { docId } = useParams<{ docId: string }>();
  const navigate = useNavigate();
  const [selectedItemsIds, setSelectedItemIds] = useState<string[]>(
    docId ? [`document_${docId}`] : []
  );
  const [activeShowItemId, setActiveShowItemId] = useState<string | undefined>(
    docId ? `document_${docId}` : undefined
  );

  const [model] = useState(() => Model.fromJson({
    global: {
      tabEnableClose: true,
      tabEnableRename: true,
      tabEnableDrag: true,
      tabSetEnableMaximize: true,
      tabSetEnableClose: true,
      tabSetEnableDeleteWhenEmpty: true,
      rootOrientationVertical: false,
    },
    layout: {
      type: 'row',
      children: [
        {
          type: 'tabset',
          id: 'tree-tabset',
          weight: 25,
          children: [{
            type: 'tab',
            id: 'tree',
            name: 'Project',
            component: 'tree',
          }],
        },
        {
          type: 'tabset',
          id: 'details-tabset',
          weight: 75,
          enableDeleteWhenEmpty: false,
          children: [],
        },
      ],
    },
  }));

  const activeDocId = activeShowItemId ? idFromDocId(activeShowItemId) : undefined;
  const activeDoc = activeDocId
    ? project.allDocuments.find(d => d.docid === activeDocId)
    : undefined;
  const previewAvailable = activeDoc ? hasPreview(activeDoc.data) : false;

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
      || !!idFromRepoId(rawShowItemId);
    const showItemId = isSpecific ? rawShowItemId : 'config';

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
    const iconSx = { fontSize: 16 };
    const icon = {
      [TabKind.Notebook]: <Code sx={iconSx} />,
      [TabKind.Document]: <Description sx={iconSx} />,
      [TabKind.Agent]: <Science sx={iconSx} />,
      [TabKind.ProjectConfig]: <Settings sx={iconSx} />,
      [TabKind.Repository]: <Source sx={iconSx} />,
      [TabKind.CentralRepository]: <FolderOpen sx={iconSx} />,
    }[kind!];
    if (icon) renderValues.leading = icon;
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
      <Layout
        model={model}
        factory={factory}
        onRenderTab={onRenderTab}
        onAction={handleAction}
      />
    </Box>
  );
};
