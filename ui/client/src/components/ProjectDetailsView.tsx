import { Paper, Typography } from '@mui/material';
import type { Project, ProjectDocument, ProjectEntire } from '../shared/types';

export const ProjectDetailsView = ({ project }: { project: ProjectEntire | null }) => {
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
      {/* <Typography>ID: {project.name}</Typography> */}
      <Typography>Documents: {project.documents?.length ?? 'N/A'}</Typography>
      {project.documents?.map(d => (
        <Typography key={d.desc.docid}>
          Name: {d.desc.datasourceName ?? 'N/A'}
        </Typography>
      ))}
      {/* <Typography>Cache Documents: {project.documents?.cache ?? 'N/A'}</Typography>
      <Typography>Toolkits: {project.toolkitCount ?? 'N/A'}</Typography> */}
    </Paper>
  );
};
