import { Box, Stack, Typography } from '@mui/material';
import { MouseEvent, ReactNode } from 'react';
import { TreeItem } from '@mui/x-tree-view';
import { useTreeViewContext, UseTreeViewExpansionSignature } from '@mui/x-tree-view/internals';
import { RenameField } from '../../elements/RenameField';
import { DetailsViewItemValue } from './DetailsViewItemValue';
import { DetailsViewItemBranchActions } from './DetailsViewItemBranchActions';
import { DeleteFieldButton } from './DeleteFieldButton';
import { ItemTypeSelector, calcItemType, ItemTypesEnum } from './ItemTypeSelector';
import { EmptyBranchLabel } from './EmptyBranchLabel';
import { DATA_FORMAT_FIELD, DESC_FIELD } from '../../shared/constants';
import { FieldDef } from './fieldDef';
import { DetailsViewItemsInArray } from './DetailsViewItemsInArray';
import { DetailsViewItemsInObject } from './DetailsViewItemsInObject';

export const keyForDetailsViewItem = (itemKey: string, parentKey?: string) => {
  return parentKey ? `${parentKey}/${itemKey}` : itemKey;
};

export const DetailsViewItem = ({
  itemKey,
  itemValue,
  setItemValue,
  setItemKey = undefined,
  isListIndex = false,
  parentKey,
  def = undefined,
  renderBeforeName = undefined,
  onRowContextMenu = undefined,
  onValueCaret = undefined,
}: {
  itemKey: string,
  itemValue: any,
  setItemValue: (newVal: any) => void,
  setItemKey?: (newKey: string | undefined) => void | undefined,
  // True for a list element, whose index is narrow, greyed and not renameable.
  isListIndex?: boolean,
  parentKey?: string,
  // This field's definition: `required` for the editor, `children` for sub-fields.
  def?: FieldDef,
  // Optional extra content before the field name.
  renderBeforeName?: (itemKey: string, parentKey: string | undefined, def?: FieldDef) => ReactNode,
  // Optional right-click handler for a row.
  onRowContextMenu?: (itemKey: string, parentKey: string | undefined, event: MouseEvent<HTMLElement>) => void,
  // Optional report of a leaf value and its caret position, for autocomplete.
  onValueCaret?: (itemKey: string, parentKey: string | undefined, value: string, caret: number | null, el: HTMLInputElement) => void,
}) => {
  const key = keyForDetailsViewItem(itemKey, parentKey);
  const itemType = calcItemType(itemValue);
  const isArray = itemType === ItemTypesEnum.array;
  const isTree = isArray || itemType === ItemTypesEnum.object;
  const level = parentKey?.split('/').length || 0;
  // Only a required leaf needs room below for its floating "required" text.
  let marginBottom = 3;
  if (!isTree && !!def?.required) {
    marginBottom = 14;
  }
  const { publicAPI } = useTreeViewContext<[UseTreeViewExpansionSignature]>();

  return (
    <TreeItem
      key={key}
      itemId={key}
      // MUI clips the label, which would cut the delete button hanging off the left.
      sx={{ '& > .MuiTreeItem-content > .MuiTreeItem-label': { overflow: 'visible' } }}
      label={(
        <Stack
          direction='row'
          spacing={1}
          justifyItems={'stretch'}
          alignItems={'center'}
          // Room below every row, so "required" text shows without moving anything.
          style={{ marginTop: 2, marginBottom }}
          sx={{
            '& .field-delete, & .field-json': { display: 'none' },
            '&:hover .field-delete, &:hover .field-json': { display: 'flex' },
          }}
          onContextMenu={event => onRowContextMenu?.(itemKey, parentKey, event)}
        >

          {renderBeforeName?.(itemKey, parentKey, def)}

          {/* The delete button sits on the name's top-left corner, over no text. */}
          <Box sx={{ position: 'relative', display: 'flex', minWidth: 0 }}>
            {/* A list index shows as [0], narrow and grey, and never renames. */}
            {isListIndex && (
              <Typography sx={{ fontFamily: 'monospace', color: 'text.secondary', whiteSpace: 'nowrap', flexShrink: 0 }}>
                {`[${itemKey}]`}
              </Typography>
            )}
            {!isListIndex && (
              <RenameField
                value={itemKey}
                setValue={setItemKey}
                labelMinWidth="100px"
                // The top-level `desc` field isn't renameable, so show a friendlier label.
                valueForView={(
                  itemKey === DESC_FIELD && !parentKey
                    ? (
                      <Typography sx={{ whiteSpace: 'nowrap', minWidth: '100px', flexShrink: 0 }}>
                        Description (desc)
                      </Typography>
                    )
                    : undefined
                )}
              />
            )}
            {setItemKey && (
              <Box
                className="field-delete"
                sx={{
                  position: 'absolute',
                  left: '-8px',
                  top: 0,
                  transform: 'translateY(-40%)',
                  zIndex: 2,
                  '& .MuiIconButton-root': { padding: '2px' },
                  '& .MuiSvgIcon-root': { fontSize: '0.8rem' },
                }}
              >
                <DeleteFieldButton
                  itemKey={itemKey}
                  setItemKey={setItemKey}
                />
              </Box>
            )}
          </Box>

          {/* The type chip picks string/number/null/object for every field,
              except dataFormat (own dropdown) and desc (hidden fields make a
              type switch unsafe). */}
          {itemKey !== DATA_FORMAT_FIELD && itemKey !== DESC_FIELD && (
            <ItemTypeSelector
              itemValue={itemValue}
              setItemValue={newVal => {
                setItemValue(newVal);
                // Switching to an object opens its new substructure.
                const newType = calcItemType(newVal);
                if (newType === ItemTypesEnum.object || newType === ItemTypesEnum.array) {
                  publicAPI.setItemExpansion({ itemId: key, shouldBeExpanded: true });
                }
              }}
            />
          )}

          {isTree && (
            <DetailsViewItemBranchActions
              itemValue={itemValue}
              setItemValue={setItemValue}
            />
          )}

          {!isTree && (
            <DetailsViewItemValue
              itemKey={itemKey}
              itemValue={itemValue}
              setItemValue={setItemValue}
              def={def}
              onCaret={onValueCaret ? (value, caret, el) => onValueCaret(itemKey, parentKey, value, caret, el) : undefined}
            />
          )}

        </Stack>
      )}
    >
      {isTree && Object.keys(itemValue).length === 0 && (
        <EmptyBranchLabel level={level} />
      )}
      {isArray && (
        <DetailsViewItemsInArray
          itemValue={itemValue}
          setItemValue={setItemValue}
          parentKey={key}
          def={def}
          renderBeforeName={renderBeforeName}
          onRowContextMenu={onRowContextMenu}
          onValueCaret={onValueCaret}
        />
      )}
      {isTree && !isArray && (
        <DetailsViewItemsInObject
          itemValue={itemValue}
          setItemValue={setItemValue}
          parentKey={key}
          def={def}
          renderBeforeName={renderBeforeName}
          onRowContextMenu={onRowContextMenu}
          onValueCaret={onValueCaret}
          isDescRoot={parentKey === undefined && itemKey === DESC_FIELD}
        />
      )}
    </TreeItem>
  )
}
