import { MoreVert } from '@mui/icons-material';
import { Divider, Popover, Stack, Typography } from '@mui/material';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { AddProjectButton } from './AddProjectButton';
import { DeleteProjectButton } from './DeleteProjectButton';
import { ExportProjectButton } from './ExportProjectButton';
import { LoadProjectButton } from './LoadProjectButton';

// Three-dots button next to the project chooser: opens a popover holding the project actions.
export const ProjectMenuButton = () => {
  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(null);

  const close = () => {
    setAnchorEl(null);
  }

  return (
    <>
      <ButtonTooltip
        title="Project actions"
        onClick={(e) => setAnchorEl(e.currentTarget)}
      >
        <MoreVert />
      </ButtonTooltip>
      <Popover
        open={Boolean(anchorEl)}
        anchorEl={anchorEl}
        onClose={close}
        anchorOrigin={{ vertical: 'bottom', horizontal: 'left' }}
      >
        <Typography variant="h6" sx={{ px: 1.5, pt: 1, pb: 0.5, color: 'text.secondary' }}>
          Project
        </Typography>
        <Divider />
        <Stack alignItems="flex-start" sx={{ px: 0.5, py: 0.5 }}>
          <AddProjectButton
            onAdded={close}
          />
          <LoadProjectButton
            onLoaded={close}
          />
          <ExportProjectButton
            onExported={close}
          />
          <DeleteProjectButton
            onDeleted={close}
          />
        </Stack>
      </Popover>
    </>
  );
};
