import { Box, Paper, Typography, useTheme } from '@mui/material';
import { useState } from 'react';
import { useParams } from 'react-router-dom';
import { DashboardHeader } from './components/header/DashboardHeader';
import { ProjectLayout } from './components/layout/ProjectLayout';
import { FetchProjects } from './io/FetchProjects';
import { ProjectAutoReload } from './io/ProjectAutoReload';
import { useProjectStore } from './stores/useProjectStore';
import { ServerConstantReader } from './stores/useServerConstants';
import { buildTabKindCss } from './shared/tabKindConfig';

export const Dashboard = () => {
  const { projectName } = useParams<{ projectName: string }>();
  const { getProject } = useProjectStore();
  const dark = useTheme().palette.mode === 'dark';
  const [treeCollapsed, setTreeCollapsed] = useState(false);
  const [layoutResetSignal, setLayoutResetSignal] = useState(0);

  const project = getProject();

  return (<>
    <ServerConstantReader />
    <FetchProjects urlProjectName={projectName} />
    <ProjectAutoReload />
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
        onResetLayout={() => setLayoutResetSignal(s => s + 1)}
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
            <Box sx={{ position: 'relative', flex: 1, height: '100%' }}>
              <style>{buildTabKindCss(dark)}</style>
              <ProjectLayout
                project={project}
                treeCollapsed={treeCollapsed}
                resetSignal={layoutResetSignal}
              />
            </Box>
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
