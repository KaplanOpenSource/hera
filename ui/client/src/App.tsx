import { Container, Stack } from '@mui/material';
import { FetcherProjectNames } from './io/fetchProject';
import { Projects } from './pages/Projects';
import { PageTitle } from './components/PageTitle';
import { CommandExecutor } from './components/CommandExecutor';
import { ProjectChooser } from './components/ProjectChooser';

export default function App() {
  return (<>
    <FetcherProjectNames />
    <Container maxWidth="lg" sx={{ py: 4 }}>
      <Stack direction={'row'} spacing={2}>
        <PageTitle />
        <ProjectChooser />
      </Stack>
      <Projects />
      <CommandExecutor />
    </Container>
  </>)
}
