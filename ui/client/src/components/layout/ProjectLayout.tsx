import { Box, Paper } from '@mui/material';
import { Actions, DockLocation, Layout, Model, TabNode } from 'flexlayout-react';
import 'flexlayout-react/style/light.css';
import { useCallback, useEffect, useState } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
import { ProjectObj } from '../../objects/ProjectObj';
import { idFromDocId } from '../../shared/idDocId';
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

  const selectedDocId = idFromDocId(selectedItemsIds[0]);
  const selectedDoc = selectedDocId
    ? project.allDocuments.find(d => d.docid === selectedDocId)
    : undefined;
  const previewAvailable = selectedDoc ? hasPreview(selectedDoc.data) : false;

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
    const tabsToRemove: string[] = [];
    model.visitNodes((node) => {
      if (node.getType() === 'tab') {
        const id = node.getId();
        if (id.startsWith('details:') || id.startsWith('preview:')) {
          tabsToRemove.push(id);
        }
      }
    });
    tabsToRemove.forEach(id => model.doAction(Actions.deleteTab(id)));

    const showItemId = selectedItemsIds[0];
    if (showItemId) {
      model.doAction(Actions.addTab(
        {
          type: 'tab',
          id: `details:${showItemId}`,
          name: detailsTabName(showItemId, project),
          component: 'details',
          config: { showItemId },
        },
        'details-tabset',
        DockLocation.CENTER,
        -1,
      ));

      if (previewAvailable && selectedDocId) {
        model.doAction(Actions.addTab(
          {
            type: 'tab',
            id: `preview:${selectedDocId}`,
            name: `Preview: ${selectedDoc!.name}`,
            component: 'preview',
            config: { docid: selectedDocId },
          },
          'details-tabset',
          DockLocation.BOTTOM,
          -1,
        ));
      }
    }
  }, [selectedItemsIds[0], previewAvailable, selectedDocId, model]);

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
      />
    </Box>
  );
};
