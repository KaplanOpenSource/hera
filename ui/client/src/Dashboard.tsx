import { Box, Paper, Typography } from '@mui/material';
import { useCallback, useEffect, useState } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
import { DetailsViewPanel } from './components/details/DetailsViewPanel';
import { DashboardHeader } from './components/header/DashboardHeader';
import { ProjectTreeView } from './components/project/ProjectTreeView';
import { SplitWithSidebar } from './elements/SplitWithSidebar';
import { FetchProjects } from './io/FetchProjects';
import { useProjectStore } from './stores/useProjectStore';
import { ServerConstantReader } from './stores/useServerConstants';

export const Dashboard = () => {
  const { projectName, docId } = useParams<{ projectName: string; docId: string }>();
  const { getProject } = useProjectStore();
  const navigate = useNavigate();
  const [selectedItemsIds, setSelectedItemIds] = useState<string[]>(
    docId ? [`document_${docId}`] : []
  );
  const [treeCollapsed, setTreeCollapsed] = useState(false);

  const project = getProject();

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

  return (<>
    <ServerConstantReader />
    <FetchProjects urlProjectName={projectName} />
    <Box
      sx={{
        height: '100vh',
        display: 'flex',
        flexDirection: 'column',
      }}
    >
      <DashboardHeader
        treeCollapsed={treeCollapsed}
        setTreeCollapsed={setTreeCollapsed}
      />

      <Box
        sx={{
          flex: 1,
          display: 'flex',
          gap: 1,
          minHeight: 0,
        }}
      >
        {project
          ? (
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
          )
          : (
            <Paper sx={{ p: 2, height: '100%', overflow: 'auto', flex: 1, minWidth: 0 }}>
              <Typography>
                No project loaded
              </Typography>
            </Paper>
          )}
      </Box>
    </Box>
  </>)
}
