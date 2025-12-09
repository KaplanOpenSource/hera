import { Alert, Box, Stack, Typography } from '@mui/material';
import { useState } from 'react';
import { CommandExecutor } from './components/CommandExecutor';
import { DetailsViewPanel } from './components/details/DetailsViewPanel';
import { PageTitle } from './components/PageTitle';
import { ProjectChooser } from './components/ProjectChooser';
import { ProjectTreeView } from './components/ProjectTreeView';
import { FetchProjects } from './io/FetchProjects';
import { useProjectStore } from './stores/useProjectStore';

export default function App() {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | undefined>(undefined);
  const { currProject } = useProjectStore();
  const [selectedItemsIds, setSelectedItemIds] = useState<string[]>([]);

  // TODO:
  // 1. when clicking on document - show it instead of project
  // 2. button to delete project - aka all its docs

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
          {currProject
            ? (
              <ProjectTreeView
                project={currProject}
                setSelectedItemIds={setSelectedItemIds}
              />
            )
            : (
              <Typography>
                No project loaded
              </Typography>
            )}
        </Box>
        <Box sx={{ width: '50%' }}>
          {currProject
            ? (
              <DetailsViewPanel
                project={currProject}
                selectedItemsIds={selectedItemsIds}
              />
            )
            : (
              <Typography>
                Select a project to see more details
              </Typography>
            )}
        </Box>
      </Box>
      <CommandExecutor />
    </Stack>
  </>)
}
