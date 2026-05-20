import { Box, Paper, Typography } from '@mui/material';
import { useState } from 'react';
import { useParams } from 'react-router-dom';
import { DashboardHeader } from './components/header/DashboardHeader';
import { ProjectLayout } from './components/layout/ProjectLayout';
import { FetchProjects } from './io/FetchProjects';
import { useProjectStore } from './stores/useProjectStore';
import { ServerConstantReader } from './stores/useServerConstants';

export const Dashboard = () => {
  const { projectName } = useParams<{ projectName: string }>();
  const { getProject } = useProjectStore();
  const [treeCollapsed, setTreeCollapsed] = useState(false);

  const project = getProject();

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
          ? <ProjectLayout project={project} treeCollapsed={treeCollapsed} />
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
