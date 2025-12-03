import { Box, Paper, Stack, Typography } from '@mui/material';
import type { Project, ProjectDocument, ProjectEntire } from '../shared/types';
import { Article, ReceiptLong } from '@mui/icons-material';

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
        <Stack key={d._id.$oid} direction='row'>
          {d.desc.datasourceName
            ? (<>
              <Article />
              <Typography>
                {d.desc.datasourceName}
              </Typography>
            </>)
            : (<>
              <ReceiptLong />
              <Typography>
                {d.type ?? 'N/A'}
              </Typography>
            </>)}
        </Stack>
      ))}
      {/* <Typography>Cache Documents: {project.documents?.cache ?? 'N/A'}</Typography>
      <Typography>Toolkits: {project.toolkitCount ?? 'N/A'}</Typography> */}
    </Paper>
  );
};
