import { Close } from '@mui/icons-material';
import { IconButton } from '@mui/material';
import { BaseEdge, EdgeLabelRenderer, EdgeProps, getBezierPath } from '@xyflow/react';
import { useState } from 'react';

export interface WorkflowRequiresEdgeData {
  hovered?: boolean;
  onRemove: () => void;
  [key: string]: unknown;
}

// Custom requires edge: shows a delete button at its midpoint while the edge (or
// the button) is hovered.
export const WorkflowRequiresEdge = ({
  id,
  sourceX,
  sourceY,
  sourcePosition,
  targetX,
  targetY,
  targetPosition,
  markerEnd,
  style,
  data,
}: EdgeProps) => {
  const [labelHover, setLabelHover] = useState(false);
  const { hovered, onRemove } = (data ?? {}) as WorkflowRequiresEdgeData;
  const [edgePath, labelX, labelY] = getBezierPath({
    sourceX, sourceY, sourcePosition, targetX, targetY, targetPosition,
  });

  return (
    <>
      <BaseEdge id={id} path={edgePath} markerEnd={markerEnd} style={style} />
      {(hovered || labelHover) && onRemove && (
        <EdgeLabelRenderer>
          <div
            className="nodrag nopan"
            onMouseEnter={() => setLabelHover(true)}
            onMouseLeave={() => setLabelHover(false)}
            style={{
              position: 'absolute',
              transform: `translate(-50%, -50%) translate(${labelX}px, ${labelY}px)`,
              pointerEvents: 'all',
              // Nodes paint after (on top of) the edge-label layer; a dataflow
              // edge's midpoint lands over the target node, so lift the button
              // above the nodes or the node swallows the click.
              zIndex: 1001,
            }}
          >
            <IconButton
              size="small"
              onClick={onRemove}
              sx={{ p: '2px', bgcolor: 'background.paper', boxShadow: 1 }}
            >
              <Close sx={{ fontSize: 14 }} />
            </IconButton>
          </div>
        </EdgeLabelRenderer>
      )}
    </>
  );
};
