import { Container } from '@mui/material';
import { FetcherProjectNames } from './io/fetchProject';
import { Projects } from './pages/Projects';
import { PageTitle } from './components/PageTitle';
import { CommandExecutor } from './components/CommandExecutor';

export default function App() {
  return (<>
    <FetcherProjectNames />
    <Container maxWidth="lg" sx={{ py: 4 }}>
      <PageTitle />
      <Projects />
      <CommandExecutor />
    </Container>
  </>)
}
