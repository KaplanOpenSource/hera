import { FormatListBulleted } from '@mui/icons-material';
import { Box, List, ListItemButton, ListItemText, Paper, Typography } from '@mui/material';
import { Panel } from '@xyflow/react';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';

// The button floats over the graph, so it carries its own surface.
const BUTTON_SX = { bgcolor: 'background.paper', boxShadow: 1, p: 0.25 };

// The canvas's list of all nodes, top-left, collapsed to one icon until opened.
// Clicking a row centres the canvas on that node and opens it for editing.
export const WorkflowNodeListPanel = ({
  nodeNames,
  selectedNode,
  onPickNode,
}: {
  nodeNames: string[],
  selectedNode?: string,
  onPickNode: (name: string) => void,
}) => {
  const [open, setOpen] = useState(false);

  let title = 'Show the node list';
  if (open) {
    title = 'Hide the node list';
  }

  return (
    <Panel position="top-left" className="nodrag">
      <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1, alignItems: 'flex-start' }}>
        <ButtonTooltip title={title} aria-label="node list" onClick={() => setOpen(!open)} sx={BUTTON_SX}>
          <FormatListBulleted sx={{ fontSize: 16 }} />
        </ButtonTooltip>
        {open && (
          <Paper sx={{ minWidth: 140, maxWidth: 240, maxHeight: 260, overflowY: 'auto' }}>
            {nodeNames.length === 0 && (
              <Typography variant="body2" color="text.secondary" sx={{ p: 1 }}>No nodes.</Typography>
            )}
            <List dense disablePadding>
              {nodeNames.map(name => (
                <ListItemButton
                  key={name}
                  selected={name === selectedNode}
                  onClick={() => onPickNode(name)}
                >
                  <ListItemText primary={name} slotProps={{ primary: { variant: 'body2', noWrap: true } }} />
                </ListItemButton>
              ))}
            </List>
          </Paper>
        )}
      </Box>
    </Panel>
  );
};
