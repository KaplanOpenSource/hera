import { Alert, Box, Container } from '@mui/material';
import { useEffect, useState } from 'react';
import { CommandExecutor } from '../components/CommandExecutor';
import { PageTitle } from '../components/PageTitle';
import { ProjectDetailsView } from '../components/ProjectDetailsView';
import { ProjectTreeView } from '../components/ProjectTreeView';
import { API_BASE } from '../shared/constants';
import type { ExecRequest, Project } from '../shared/types';
import { execPython } from '../io/execPython';

export const Projects = () => {
  const [projects, setProjects] = useState<Project[]>([]);
  const [selectedProject, setSelectedProject] = useState<Project | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | undefined>(undefined);

  useEffect(() => {
    const fetchProjectNames = async () => {
      setLoading(true);
      setError(undefined);
      const { data, problem } = await execPython(`
from hera.datalayer.project import getProjectList;
result = [{"id": "p-" + str(i), "name": proj} for i, proj in enumerate(getProjectList())]
          `)
      if (data) {
        setProjects(data || []);
      } else {
        setError(problem);
      }
      setLoading(false);
    };
    fetchProjectNames();
  }, []);

  const handleProjectSelect = (project: Project) => {
    setSelectedProject(project);
  };

  const fetchProjectDetails = async (projectId: string) => {
    const project = projects.find((p) => p.id === projectId);
    if (project && project.documents) {
      return; // Details already fetched
    }

    try {
      const payload: ExecRequest = {
        code: `result = next((p for p in MOCK_PROJECTS), None)`,
        // code: `result = next((p for p in MOCK_PROJECTS if p.id == "${projectId}"), None)`,
      };
      const r = await fetch(`${API_BASE}/exec`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload),
      });
      const data = await r.json();
      if (data) {
        setProjects((prevProjects) =>
          prevProjects.map((p) => (p.id === projectId ? { ...p, ...data } : p))
        );
      }
    } catch (e: any) {
      console.error('Failed to fetch project details:', e);
    }
  };

  return (
    <Container maxWidth="lg" sx={{ py: 4 }}>
      <PageTitle />
      {error && (
        <Box sx={{ mb: 2 }}>
          <Alert severity="error">{error}</Alert>
        </Box>
      )}
      <Box sx={{ display: 'flex', gap: 2 }}>
        <Box sx={{ width: '50%' }}>
          <ProjectTreeView
            projects={projects}
            onProjectSelect={handleProjectSelect}
            onProjectExpand={fetchProjectDetails}
          />
        </Box>
        <Box sx={{ width: '50%' }}>
          <ProjectDetailsView project={selectedProject} />
        </Box>
      </Box>
      <CommandExecutor />
    </Container>
  );
};


