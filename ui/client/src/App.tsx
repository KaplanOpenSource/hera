import { Alert, Box, Stack } from '@mui/material';
import { useEffect, useState } from 'react';
import { CommandExecutor } from './components/CommandExecutor';
import { PageTitle } from './components/PageTitle';
import { ProjectChooser } from './components/ProjectChooser';
import { ProjectDetailsView } from './components/ProjectDetailsView';
import { ProjectTreeView } from './components/ProjectTreeView';
import { FetchProjects } from './io/FetchProjects';
import { useProjectStore } from './stores/useProjectStore';

export default function App() {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | undefined>(undefined);
  const { currProject } = useProjectStore();

  return (<>
    <FetchProjects />
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
