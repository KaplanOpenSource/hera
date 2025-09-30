import { Paper, Typography } from '@mui/material';
import type { Project } from '../shared/types';

export const ProjectDetailsView = ({ project }: { project: Project | null }) => {
  if (!project) {
    return (
      <Paper sx={{ p: 2, height: '100%' }}>
        <Typography variant="h6">Project Details</Typography>
        <Typography>Select a project to see more details.</Typography>
      </Paper>
    );
  }

  return (
    <Paper sx={{ p: 2, height: '100%' }}>
      <Typography variant="h6">{project.name}</Typography>
      <Typography>ID: {project.id}</Typography>
      <Typography>Sim Documents: {project.documents?.sim ?? 'N/A'}</Typography>
      <Typography>Measure Documents: {project.documents?.measure ?? 'N/A'}</Typography>
      <Typography>Cache Documents: {project.documents?.cache ?? 'N/A'}</Typography>
      <Typography>Toolkits: {project.toolkitCount ?? 'N/A'}</Typography>
    </Paper>
  );
};
