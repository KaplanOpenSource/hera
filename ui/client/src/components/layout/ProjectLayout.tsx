import { Visibility } from '@mui/icons-material';
import { Box, Paper } from '@mui/material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { useCallback, useEffect, useState } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
import { SplitWithSidebar } from '../../elements/SplitWithSidebar';
import { ProjectObj } from '../../objects/ProjectObj';
import { DetailsViewPanel } from '../details/DetailsViewPanel';
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
      <Box sx={{ height: '100%', position: 'relative' }}>
        <Paper sx={{
          height: '100%',
          overflow: 'hidden',
        }}>
          <DetailsViewPanel
            project={project}
            showItemId={selectedItemsIds[0]}
            previewHidden={previewHidden}
            setPreviewHidden={setPreviewHidden}
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
    </SplitWithSidebar>
  );
};
