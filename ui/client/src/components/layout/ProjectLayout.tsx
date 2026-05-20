import { Visibility, VisibilityOff } from '@mui/icons-material';
import { Box, Paper } from '@mui/material';
import { useCallback, useEffect, useState } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
import { Panel, Group as PanelGroup, Separator as PanelResizeHandle, usePanelRef } from 'react-resizable-panels';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { SplitWithSidebar } from '../../elements/SplitWithSidebar';
import { ProjectObj } from '../../objects/ProjectObj';
import { idFromDocId } from '../../shared/idDocId';
import { DetailsViewPanel } from '../details/DetailsViewPanel';
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
  const [previewHidden, setPreviewHidden] = useState(false);
  const previewPanelRef = usePanelRef();

  const selectedDocId = idFromDocId(selectedItemsIds[0]);
  const selectedDoc = selectedDocId
    ? project.allDocuments.find(d => d.docid === selectedDocId)
    : undefined;
  const previewAvailable = selectedDoc ? hasPreview(selectedDoc.data) : false;

  useEffect(() => {
    setPreviewHidden(false);
  }, [selectedDocId]);

  useEffect(() => {
    if (previewHidden) {
      previewPanelRef.current?.collapse();
    } else {
      previewPanelRef.current?.expand();
    }
  }, [previewHidden]);

  useEffect(() => {
    if (docId && project?.documentIds.has(docId)) {
      setSelectedItemIds([`document_${docId}`]);
    } else {
      setSelectedItemIds([]);
    }
  }, [project?.name])

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

  return (
    <SplitWithSidebar
      collapsed={treeCollapsed}
      sidebar={
        <Paper sx={{ p: 2, height: '100%', overflow: 'auto' }}>
          <ProjectTreeView
            project={project}
            selectedItemsIds={selectedItemsIds}
            setSelectedItemIds={handleSetSelectedItemIds}
          />
        </Paper>
      }
    >
      {previewAvailable
        ? (
          <PanelGroup
            orientation="vertical"
            onLayoutChanged={(layout) => {
              if (layout['preview-panel'] === 0) {
                setPreviewHidden(true);
              }
            }}
          >
            <Panel defaultSize={50} minSize={20}>
              <Box sx={{ height: '100%', position: 'relative' }}>
                <Paper sx={{ height: '100%', overflow: 'hidden' }}>
                  <DetailsViewPanel
                    project={project}
                    showItemId={selectedItemsIds[0]}
                  />
                </Paper>
                {previewHidden && (
                  <Box sx={{ position: 'absolute', bottom: 4, right: 4 }}>
                    <ButtonTooltip
                      title="Show preview"
                      onClick={() => setPreviewHidden(false)}
                    >
                      <Visibility sx={{ fontSize: 14 }} />
                    </ButtonTooltip>
                  </Box>
                )}
              </Box>
            </Panel>

            <PanelResizeHandle
              style={{
                height: 4,
                cursor: 'row-resize',
                backgroundColor: '#e0e0e0',
                outline: 'none',
              }}
            />

            <Panel
              id="preview-panel"
              panelRef={previewPanelRef}
              defaultSize={50}
              minSize={5}
              collapsible
            >
              <Box sx={{ height: '100%', position: 'relative' }}>
                <Box sx={{ position: 'absolute', top: 4, right: 4, zIndex: 1000 }}>
                  <ButtonTooltip
                    title="Hide preview"
                    onClick={() => setPreviewHidden(true)}
                    sx={{
                      backgroundColor: 'white',
                      '&:hover': { backgroundColor: '#eee' },
                    }}
                  >
                    <VisibilityOff sx={{ fontSize: 14 }} />
                  </ButtonTooltip>
                </Box>
                <PreviewPanel docid={selectedDocId!} />
              </Box>
            </Panel>
          </PanelGroup>
        )
        : (
          <Box sx={{ height: '100%', position: 'relative' }}>
            <Paper sx={{ height: '100%', overflow: 'hidden' }}>
              <DetailsViewPanel
                project={project}
                showItemId={selectedItemsIds[0]}
              />
            </Paper>
            {previewHidden && (
              <Box sx={{ position: 'absolute', bottom: 4, right: 4 }}>
                <ButtonTooltip
                  title="Show preview"
                  onClick={() => setPreviewHidden(false)}
                >
                  <Visibility sx={{ fontSize: 14 }} />
                </ButtonTooltip>
              </Box>
            )}
          </Box>
        )
      }
    </SplitWithSidebar>
  );
};
