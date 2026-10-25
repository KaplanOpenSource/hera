import { MenuItem, TextField } from '@mui/material';
import { WorkflowChunk } from '../../../io/runWorkflow';
import { chunkHasContent, isTaskChunk } from './WorkflowChunkLog';

// Picks which task's output to jump to, and shows the one currently in view.
export const WorkflowNodeSelect = ({
  chunks,
  currentIndex,
  onPick,
}: {
  chunks: WorkflowChunk[],
  currentIndex?: number,
  onPick: (index: number) => void,
}) => {
  // Only task chunks that actually render a card can be scrolled to.
  const options = chunks
    .map((chunk, index) => { return { index, chunk }; })
    .filter((option) => { return isTaskChunk(option.chunk.name) && chunkHasContent(option.chunk); });

  const value = options.some(o => o.index === currentIndex) ? String(currentIndex) : '';

  return options.length > 0 && (
    <TextField
      select
      size="small"
      label="Node"
      value={value}
      onChange={(e) => { return onPick(Number(e.target.value)); }}
      sx={{ minWidth: 140, bgcolor: 'background.paper', boxShadow: 1, borderRadius: 1 }}
      slotProps={{ select: { MenuProps: { PaperProps: { sx: { maxHeight: 320 } } } } }}
    >
      {options.map((option) => {
        return <MenuItem key={option.index} value={String(option.index)}>{option.chunk.name}</MenuItem>;
      })}
    </TextField>
  );
};
