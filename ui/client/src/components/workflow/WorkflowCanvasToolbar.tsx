import { Add, AutoAwesome } from '@mui/icons-material';
import { Box, Divider, Menu, MenuItem } from '@mui/material';
import { Panel } from '@xyflow/react';
import { ReactNode, useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { WorkflowBlock } from '../../shared/types';
import { workflowTemplates } from './workflowTemplates';

// The buttons float over the graph, so each carries its own surface.
const BUTTON_SX = { bgcolor: 'background.paper', boxShadow: 1, p: 0.25 };

// The canvas's top-right controls: add a node, pick a starter template, plus
// whatever the page passes in (e.g. the run button). Owns the templates menu.
export const WorkflowCanvasToolbar = ({
  actionButtons,
  onAddNode,
  onApplyTemplate,
}: {
  actionButtons?: ReactNode,
  onAddNode: () => void,
  onApplyTemplate: (block: WorkflowBlock) => void,
}) => {
  const [templatesAnchor, setTemplatesAnchor] = useState<HTMLElement | null>(null);

  const applyTemplate = (block: WorkflowBlock) => {
    onApplyTemplate(block);
    setTemplatesAnchor(null);
  };

  return (
    <Panel position="top-right" style={{ marginRight: 24 }}>
      <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1, alignItems: 'flex-end' }}>
        <ButtonTooltip title="Add node" onClick={onAddNode} sx={BUTTON_SX}>
          <Add sx={{ fontSize: 16 }} />
        </ButtonTooltip>
        <ButtonTooltip title="Templates" onClick={e => setTemplatesAnchor(e.currentTarget)} sx={BUTTON_SX}>
          <AutoAwesome sx={{ fontSize: 16 }} />
        </ButtonTooltip>
        {actionButtons}
      </Box>
      <Menu anchorEl={templatesAnchor} open={!!templatesAnchor} onClose={() => setTemplatesAnchor(null)}>
        {workflowTemplates.map(t => (
          <MenuItem key={t.id} onClick={() => applyTemplate(t.block)}>
            {t.label}
          </MenuItem>
        ))}
        <Divider />
        <MenuItem onClick={() => setTemplatesAnchor(null)}>Cancel</MenuItem>
      </Menu>
    </Panel>
  );
};
