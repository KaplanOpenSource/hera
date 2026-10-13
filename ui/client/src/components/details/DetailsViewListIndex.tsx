import { Stack, Typography } from '@mui/material';
import { DragIndicator } from '@mui/icons-material';
import { DraggableAttributes } from '@dnd-kit/core';
import { SyntheticListenerMap } from '@dnd-kit/core/dist/hooks/utilities';

// A list element's index, next to a grip. Both are the handle for reordering.
export const DetailsViewListIndex = ({
  index,
  attributes,
  listeners,
}: {
  index: number,
  attributes: DraggableAttributes,
  listeners?: SyntheticListenerMap,
}) => {
  return (
    <Stack
      direction="row"
      alignItems="center"
      spacing={0.25}
      data-testid={`list-index-${index}`}
      onClick={e => e.stopPropagation()}
      sx={{
        color: 'text.secondary',
        flexShrink: 0,
        cursor: 'grab',
        touchAction: 'none',
      }}
      {...attributes}
      {...listeners}
    >
      <DragIndicator sx={{ fontSize: '1rem' }} />
      <Typography sx={{ fontFamily: 'monospace' }}>
        {index}
      </Typography>
    </Stack>
  );
};
