import { Alert, Box, Container, Stack } from '@mui/material';
import { FetcherProjectNames } from './io/fetchProject';
import { Projects } from './pages/Projects';
import { PageTitle } from './components/PageTitle';
import { CommandExecutor } from './components/CommandExecutor';
import { ProjectChooser } from './components/ProjectChooser';
import { useState } from 'react';
import { ProjectDetailsView } from './components/ProjectDetailsView';
import { useProjectStore } from './stores/useProjectStore';

export default function App() {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | undefined>(undefined);
  const { currProject } = useProjectStore();

  return (<>
    <FetcherProjectNames />
    <Stack spacing={2} margin={2}>
      <Stack direction={'row'} spacing={2}>
        <PageTitle />
        <ProjectChooser />
      </Stack>

      {/* $$$$ TODO:
      instead of showing at tree of Project state, show the currProject with type ProjectEntire */}

      {error && (
        <Box sx={{ mb: 2 }}>
          <Alert severity="error">{error}</Alert>
        </Box>
      )}
      <Box sx={{ display: 'flex', gap: 2, height: '80vh' }}>
        <Box sx={{ width: '50%' }}>
          {/* <ProjectTreeView
            projects={projects}
            onProjectSelect={handleProjectSelect}
            onProjectExpand={fetchProjectDetails}
          /> */}
          project tree view
        </Box>
        <Box sx={{ width: '50%' }}>
          {/* project details */}
          <ProjectDetailsView project={currProject} />
        </Box>
      </Box>
      <CommandExecutor />
    </Stack>
  </>)
}
