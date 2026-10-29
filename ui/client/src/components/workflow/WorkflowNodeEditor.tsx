import { Autocomplete, Box, InputBase, Stack, TextField, Typography } from '@mui/material';
import { useState } from 'react';
import { WorkflowNode } from '../../shared/types';
import { normalizeRequires } from '../../shared/workflow';
import { keyForDetailsViewItem } from '../details/DetailsViewItem';
import { NodeCatalogEntry, nodeOutputNames, nodeTypeGroup, nodeTypeIssue, paramsFieldDef } from './nodeCatalog';
import { paramsOnTypeChange } from './nodeTypeParams';
import { INPUT_PARAMETERS_KEY, WorkflowNodeInputs } from './WorkflowNodeInputs';
import { WorkflowNodeOutputs } from './WorkflowNodeOutputs';
import { WorkflowRequiresLine } from './WorkflowRequiresLine';
import { useNodeParamsValidation } from './useNodeParamsValidation';

// A node's editor (name, type, parameters, outputs), on the canvas or in the pop-up.
export const WorkflowNodeEditor = ({
  name,
  node,
  catalog,
  onRename,
  onChange,
  onFieldContextMenu,
  onFieldInlineEdit,
  withHandles = true,
}: {
  name: string,
  node: WorkflowNode,
  catalog: NodeCatalogEntry[],
  onRename: (newName: string) => void,
  onChange: (node: WorkflowNode) => void,
  onFieldContextMenu: (param: string, x: number, y: number, caret?: number) => void,
  onFieldInlineEdit: (param: string, value: string, caret: number | null, el: HTMLInputElement) => void,
  // False off the canvas (the pop-up), where no line can start or land.
  withHandles?: boolean,
}) => {
  const [draft, setDraft] = useState(name);

  const params = node.Execution?.input_parameters ?? {};

  // Expanded two levels by default: input_parameters and each parameter under it,
  // so nested parameter values are visible without a click.
  const [expandedItems, setExpandedItems] = useState<string[]>(() => {
    const inputsKey = keyForDetailsViewItem(INPUT_PARAMETERS_KEY);
    return [inputsKey, ...Object.keys(params).map(param => keyForDetailsViewItem(param, inputsKey))];
  });

  const typeOptions = catalog.map(entry => entry.type);
  const typeIssue = nodeTypeIssue(node, catalog);
  // Hermes's own complaint about the parameter values, shown under the node.
  const paramsIssue = useNodeParamsValidation(node);
  const paramsDef = paramsFieldDef(node, catalog);
  const outputs = nodeOutputNames(node, catalog);

  // Free-form typing keeps the type as-is (custom types stay allowed); picking a
  // known type also seeds its parameters from the catalog.
  const setType = (type: string) => onChange({ ...node, type });
  const pickType = (type: string) => {
    const entry = catalog.find(e => e.type === type);
    onChange({
      ...node,
      type,
      Execution: { ...node.Execution, input_parameters: paramsOnTypeChange(params, entry) },
    });
  };

  const commit = () => {
    const next = draft.trim();
    if (next && next !== name) {
      onRename(next);
    } else {
      setDraft(name);
    }
  };

  return (
    <>
      <InputBase
        className="nodrag"
        value={draft}
        onChange={(e) => setDraft(e.target.value)}
        onBlur={commit}
        onKeyDown={(e) => { if (e.key === 'Enter') { (e.target as HTMLInputElement).blur(); } }}
        inputProps={{ style: { padding: 0 }, 'aria-label': 'node name' }}
        sx={{ fontSize: 13, fontWeight: 600 }}
      />
      <Box className="nodrag" sx={{ mt: 1 }}>
        <Autocomplete
          className="nodrag"
          freeSolo
          size="small"
          options={typeOptions}
          groupBy={(option) => nodeTypeGroup(option)}
          inputValue={node.type ?? ''}
          onInputChange={(_e, value, reason) => {
            if (reason === 'input') {
              setType(value);
            } else if (reason === 'clear') {
              setType('');
            }
          }}
          onChange={(_e, value) => pickType(typeof value === 'string' ? value : value ?? '')}
          renderInput={(inputParams) => <TextField {...inputParams} label="type" fullWidth />}
        />
        <Box sx={{ mt: 0.5 }}>
          <WorkflowRequiresLine requires={normalizeRequires(node.requires)} />
        </Box>
        <Stack spacing={1} sx={{ mt: 1, alignItems: 'stretch' }}>
          <WorkflowNodeInputs
            nodeName={name}
            params={params}
            paramsDef={paramsDef}
            expandedItems={expandedItems}
            onExpandedItemsChange={setExpandedItems}
            onChangeParams={(newVal) => onChange({
              ...node,
              Execution: { ...node.Execution, input_parameters: newVal },
            })}
            onFieldContextMenu={onFieldContextMenu}
            onFieldInlineEdit={onFieldInlineEdit}
            withHandles={withHandles}
          />
          {outputs.length > 0 && (
            <WorkflowNodeOutputs nodeName={name} outputs={outputs} withHandles={withHandles} />
          )}
        </Stack>
        {typeIssue && (
          <Typography
            className="nodrag"
            variant="caption"
            color="warning.main"
            sx={{ display: 'block', mt: 0.5, userSelect: 'text' }}
          >
            {typeIssue}
          </Typography>
        )}
        {paramsIssue && (
          <Typography
            className="nodrag"
            variant="caption"
            color="error.main"
            sx={{ display: 'block', mt: 0.5, userSelect: 'text' }}
          >
            {paramsIssue}
          </Typography>
        )}
      </Box>
    </>
  );
};
