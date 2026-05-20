import { Paper } from '@mui/material';
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
      <Paper sx={{
        height: '100%',
        overflow: 'hidden',
      }}>
        <DetailsViewPanel
          project={project}
          showItemId={selectedItemsIds[0]}
        />
      </Paper>
    </SplitWithSidebar>
  );
};
