import { MoreVert } from '@mui/icons-material';
import { Popover, Stack } from '@mui/material';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { AddDocumentButton } from './AddDocumentButton';

// "Actions" button on the project node: opens a popover holding the project-level
// actions (currently just Add Document).
export const ProjectActionsButton = ({
  onDocumentCreated,
}: {
  onDocumentCreated?: (docOid: string) => void,
}) => {
  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(null);

  return (
    <>
      <ButtonTooltip
        title="Actions"
        onClick={(e) => setAnchorEl(e.currentTarget)}
      >
        <MoreVert />
      </ButtonTooltip>
      <Popover
        open={Boolean(anchorEl)}
        anchorEl={anchorEl}
        onClose={() => setAnchorEl(null)}
        anchorOrigin={{ vertical: 'bottom', horizontal: 'left' }}
        onClick={(e) => e.stopPropagation()}
      >
        <Stack direction="row" alignItems="center" sx={{ px: 0.5, py: 0.25 }}>
          <AddDocumentButton
            onDocumentCreated={(oid) => {
              onDocumentCreated?.(oid);
              setAnchorEl(null);
            }}
          />
        </Stack>
      </Popover>
    </>
  );
};
