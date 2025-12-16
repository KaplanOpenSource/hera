import { Alert, Box, Stack, Typography } from '@mui/material';
import { useEffect, useState } from 'react';
import { CommandExecutor } from './components/CommandExecutor';
import { DetailsViewPanel } from './components/details/DetailsViewPanel';
import { PageTitle } from './components/PageTitle';
import { ProjectChooser } from './components/ProjectChooser';
import { ProjectTreeView } from './components/ProjectTreeView';
import { FetchProjects } from './io/FetchProjects';
import { useProjectStore } from './stores/useProjectStore';
import { ServerConstantReader } from './stores/useServerConstants';

export default function App() {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | undefined>(undefined);
  const { currProject } = useProjectStore();
  const [selectedItemsIds, setSelectedItemIds] = useState<string[]>([]);

  useEffect(() => {
    setSelectedItemIds([]);
  }, [currProject?.name])

  return (<>
    <ServerConstantReader />
    <FetchProjects />
    <Stack spacing={2} margin={2}>
      <Stack direction={'row'}>
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
