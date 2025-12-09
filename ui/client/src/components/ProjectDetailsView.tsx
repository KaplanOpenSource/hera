import { Article, ReceiptLong } from '@mui/icons-material';
import { Paper, Stack, Typography } from '@mui/material';
import type { ProjectEntire } from '../shared/types';
import { useEffect } from 'react';
import { idFromDocId } from '../shared/idDocId';

export const ProjectDetailsView = ({
  project,
  selectedItemsIds,
}: {
  project: ProjectEntire | null,
  selectedItemsIds: string[],
}) => {
  if (!project) {
    return (
      <Paper sx={{ p: 2, height: '100%' }}>
        <Typography variant="h6">Project Details</Typography>
        <Typography>Select a project to see more details.</Typography>
      </Paper>
    );
  }

  useEffect(() => {
    console.log(idFromDocId(selectedItemsIds[0]));
  }, [selectedItemsIds])

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
