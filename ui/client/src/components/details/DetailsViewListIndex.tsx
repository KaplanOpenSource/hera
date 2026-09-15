import { useState } from 'react';
import { Typography, useTheme } from '@mui/material';

// A list element's index. It is the drag handle for reordering the list.
export const DetailsViewListIndex = ({
  index,
  onMove,
}: {
  index: number,
  onMove: (from: number, to: number) => void,
}) => {
  const [isDropTarget, setIsDropTarget] = useState(false);
  const theme = useTheme();

  // Drag a copy of the whole row, shown as a raised paper.
  const setDragImage = (e: React.DragEvent) => {
    const row = (e.target as HTMLElement).closest('.MuiTreeItem-content');
    if (!row || !e.dataTransfer.setDragImage) {
      return;
    }
    const paper = document.createElement('div');
    paper.appendChild(row.cloneNode(true));
    paper.style.position = 'fixed';
    paper.style.top = '-1000px';
    paper.style.width = `${row.clientWidth}px`;
    paper.style.padding = '4px 8px';
    paper.style.boxSizing = 'border-box';
    paper.style.background = theme.palette.background.paper;
    paper.style.borderRadius = `${theme.shape.borderRadius}px`;
    paper.style.boxShadow = theme.shadows[8];
    document.body.appendChild(paper);
    e.dataTransfer.setDragImage(paper, 16, 16);
    setTimeout(() => paper.remove(), 0);
  };

  const onDrop = (e: React.DragEvent) => {
    e.preventDefault();
    e.stopPropagation();
    setIsDropTarget(false);
    const from = Number(e.dataTransfer.getData('text/plain'));
    if (Number.isInteger(from)) {
      onMove(from, index);
    }
  };

  // The line marks where the dragged element will land.
  const sx: any = {
    fontFamily: 'monospace',
    color: 'text.secondary',
    whiteSpace: 'nowrap',
    flexShrink: 0,
    cursor: 'grab',
    borderTop: '2px solid transparent',
  };
  if (isDropTarget) {
    sx.borderTop = '2px solid';
    sx.borderColor = 'primary.main';
  }

  return (
    <Typography
      draggable
      data-testid={`list-index-${index}`}
      onDragStart={e => {
        e.dataTransfer.setData('text/plain', String(index));
        e.dataTransfer.effectAllowed = 'move';
        setDragImage(e);
      }}
      onDragOver={e => {
        e.preventDefault();
        e.dataTransfer.dropEffect = 'move';
        setIsDropTarget(true);
      }}
      onDragLeave={() => setIsDropTarget(false)}
      onDragEnd={() => setIsDropTarget(false)}
      onDrop={onDrop}
      onClick={e => e.stopPropagation()}
      sx={sx}
    >
      {`[${index}]`}
    </Typography>
  );
};
