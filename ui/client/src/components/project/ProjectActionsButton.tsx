import { MoreVert } from '@mui/icons-material';
import { Divider, Popover, Stack, Typography } from '@mui/material';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { CopyProjectButton } from '../header/CopyProjectButton';
import { DeleteProjectButton } from '../header/DeleteProjectButton';
import { AddDocumentButton } from './AddDocumentButton';
import { DeleteSelectedButton } from './DeleteSelectedButton';
import { DetectNotebooksButton } from './DetectNotebooksButton';

// "Actions" button on the project node: opens a popover holding the project-level actions.
export const ProjectActionsButton = ({
  selectedIds = [],
  onSelectDocument,
}: {
  selectedIds?: string[],
  // Sets the tree selection to a document (after creating one) or clears it (after deleting).
  onSelectDocument?: (docOid?: string) => void,
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
        <Typography variant="h6" sx={{ px: 1.5, pt: 1, pb: 0.5, color: 'text.secondary' }}>
          Actions
        </Typography>
        <Divider />
        <Stack alignItems="flex-start" sx={{ px: 0.5, py: 0.5 }}>
          <AddDocumentButton
            onDocumentCreated={(oid) => {
              onSelectDocument?.(oid);
              setAnchorEl(null);
            }}
          />
          <DetectNotebooksButton
            onDetected={() => setAnchorEl(null)}
          />
          <DeleteSelectedButton
            selectedIds={selectedIds}
            onDeleted={() => {
              onSelectDocument?.(undefined);
              setAnchorEl(null);
            }}
          />
          <CopyProjectButton
            onCopied={() => setAnchorEl(null)}
          />
          <DeleteProjectButton
            onDeleted={() => setAnchorEl(null)}
          />
        </Stack>
      </Popover>
    </>
  );
};
