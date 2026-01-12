import { Alert, Box, Stack, Typography } from '@mui/material';
import { useEffect, useState } from 'react';
import { CommandExecutor } from './components/CommandExecutor';
import { DetailsViewPanel } from './components/details/DetailsViewPanel';
import { PageTitle } from './components/header/PageTitle';
import { ProjectChooser } from './components/header/ProjectChooser';
import { ProjectTreeView } from './components/project/ProjectTreeView';
import { FetchProjects } from './io/FetchProjects';
import { useProjectStore } from './stores/useProjectStore';
import { ServerConstantReader } from './stores/useServerConstants';
import { CommitIdShower } from './components/header/CommitIdShower';

export default function App() {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | undefined>(undefined);
  const { getProject } = useProjectStore();
  const [selectedItemsIds, setSelectedItemIds] = useState<string[]>([]);

  const project = getProject();

  useEffect(() => {
    setSelectedItemIds([]);
  }, [project?.name])

  return (<>
    <ServerConstantReader />
    <FetchProjects />
    <Stack spacing={2} margin={2}>
      <Stack direction={'row'}>
        <PageTitle />
        <ProjectChooser />
        <CommitIdShower />
      </Stack>
      {error && (
        <Box sx={{ mb: 2 }}>
          <Alert severity="error">{error}</Alert>
        </Box>
      )}
      <Box sx={{ display: 'flex', gap: 2, height: '80vh' }}>
        <Box sx={{ width: '50%' }}>
          {project
            ? (
              <ProjectTreeView
                project={project}
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
          {project
            ? (
              <DetailsViewPanel
                project={project}
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
