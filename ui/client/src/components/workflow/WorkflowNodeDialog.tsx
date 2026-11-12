import { Close } from '@mui/icons-material';
import { Dialog, DialogContent, DialogTitle, IconButton } from '@mui/material';
import { WorkflowNode } from '../../shared/types';
import { NodeCatalogEntry } from './nodeCatalog';
import { WorkflowNodeEditor } from './WorkflowNodeEditor';

// A workflow node's editor in a pop-up, away from the canvas, so it has no line handles.
export const WorkflowNodeDialog = ({
  name,
  node,
  catalog,
  onClose,
  onRename,
  onChange,
  onFieldContextMenu,
  onFieldInlineEdit,
}: {
  name: string,
  node: WorkflowNode,
  catalog: NodeCatalogEntry[],
  onClose: () => void,
  onRename: (newName: string) => void,
  onChange: (node: WorkflowNode) => void,
  onFieldContextMenu: (param: string, x: number, y: number, caret?: number) => void,
  onFieldInlineEdit: (param: string, value: string, caret: number | null, el: HTMLInputElement) => void,
}) => {
  return (
    <Dialog open onClose={onClose} maxWidth="md" fullWidth>
      <DialogTitle sx={{ pr: 6 }}>
        {name}
        <IconButton aria-label="close" onClick={onClose} sx={{ position: 'absolute', top: 8, right: 8 }}>
          <Close />
        </IconButton>
      </DialogTitle>
      <DialogContent>
        <WorkflowNodeEditor
          name={name}
          node={node}
          catalog={catalog}
          onRename={onRename}
          onChange={onChange}
          onFieldContextMenu={onFieldContextMenu}
          onFieldInlineEdit={onFieldInlineEdit}
          withHandles={false}
        />
      </DialogContent>
    </Dialog>
  );
};
