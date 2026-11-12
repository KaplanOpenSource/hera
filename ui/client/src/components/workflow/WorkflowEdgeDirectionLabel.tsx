import { orange } from '@mui/material/colors';

// One end label: says how the node at that end relates to the node at the other.
export const WorkflowEdgeDirectionLabel = ({
  text,
  x,
  y,
}: {
  text: string,
  x: number,
  y: number,
}) => {
  return (
    <div
      className="nodrag nopan"
      style={{
        position: 'absolute',
        transform: `translate(-50%, -50%) translate(${x}px, ${y}px)`,
        fontSize: 9,
        lineHeight: 1.2,
        padding: '0 3px',
        borderRadius: 3,
        whiteSpace: 'nowrap',
        background: orange[700],
        color: 'white',
        zIndex: 1001,
      }}
    >
      {text}
    </div>
  );
};
