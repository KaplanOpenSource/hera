import { Alert, Box, Container, Stack } from '@mui/material';
import { FetcherProjectNames } from './io/fetchProject';
import { Projects } from './pages/Projects';
import { PageTitle } from './components/PageTitle';
import { CommandExecutor } from './components/CommandExecutor';
import { ProjectChooser } from './components/ProjectChooser';
import { useState } from 'react';
import { ProjectDetailsView } from './components/ProjectDetailsView';
import { useProjectStore } from './stores/useProjectStore';
import { ProjectTreeView } from './components/ProjectTreeView';

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
      {error && (
        <Box sx={{ mb: 2 }}>
          <Alert severity="error">{error}</Alert>
        </Box>
      )}
      <Box sx={{ display: 'flex', gap: 2, height: '80vh' }}>
        <Box sx={{ width: '50%' }}>
          <ProjectTreeView
            project={currProject}
          />
        </Box>
        <Box sx={{ width: '50%' }}>
          <ProjectDetailsView
            project={currProject}
          />
        </Box>
      </Box>
      <CommandExecutor />
    </Stack>
  </>)
}
