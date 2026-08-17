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
import { DATA_FORMAT_FIELD, DESC_FIELD, FILES_DIRECTORY_FIELD } from '../../shared/constants';
import { FieldDef } from './fieldDef';

export const keyForDetailsViewItem = (itemKey: string, parentKey?: string) => {
  return parentKey ? `${parentKey}/${itemKey}` : itemKey;
};

export const DetailsViewItem = ({
  itemKey,
  itemValue,
  setItemValue,
  setItemKey = undefined,
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
  parentKey?: string,
  // Definition of this field: `required` for the value editor, `children` for
  // the sub-fields below it.
  def?: FieldDef,
  // Optional slot rendered before the field name, given the row's key info and
  // its def (e.g. the workflow editor renders a source dot wrapped in a connection
  // handle). Nothing is rendered when not provided. Passed down the tree.
  renderBeforeName?: (itemKey: string, parentKey: string | undefined, def?: FieldDef) => ReactNode,
  // Optional right-click handler for a row, given the row's key info (e.g. the
  // workflow editor opens a per-field menu). Passed down the tree.
  onRowContextMenu?: (itemKey: string, parentKey: string | undefined, event: MouseEvent<HTMLElement>) => void,
  // Optional caret/value reporter for a row's leaf value editor, tagged with the
  // row's key info (e.g. the workflow editor drives inline autocomplete). Passed
  // down the tree.
  onValueCaret?: (itemKey: string, parentKey: string | undefined, value: string, caret: number | null, el: HTMLInputElement) => void,
}) => {
  const key = keyForDetailsViewItem(itemKey, parentKey);
  const isTree = calcItemType(itemValue) === ItemTypesEnum.object;
  const level = parentKey?.split('/').length || 0;
  // Only required leaf fields need bottom room for their floating "required"
  // helper text; every other row stays compact.
  let marginBottom = 3;
  if (!isTree && !!def?.required) {
    marginBottom = 14;
  }
  const { publicAPI } = useTreeViewContext<[UseTreeViewExpansionSignature]>();

  return (
    <TreeItem
      key={key}
      itemId={key}
      label={(
        <Stack
          direction='row'
          spacing={1}
          justifyItems={'stretch'}
          alignItems={'center'}
          // Bottom space on every row reserves room for a field's "required"
          // helper text, so it shows without moving anything.
          style={{ marginTop: 2, marginBottom }}
          sx={{
            '& .field-delete, & .field-json': { display: 'none' },
            '&:hover .field-delete, &:hover .field-json': { display: 'flex' },
          }}
          onContextMenu={event => onRowContextMenu?.(itemKey, parentKey, event)}
        >

          {renderBeforeName?.(itemKey, parentKey, def)}

          {/* The delete button overlays the right side of the name, so showing
              it on hover never changes the row's size. */}
          <Box sx={{ position: 'relative', display: 'flex', minWidth: 0 }}>
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
            {setItemKey && (
              <Box
                className="field-delete"
                sx={{
                  position: 'absolute',
                  right: 0,
                  top: '50%',
                  transform: 'translateY(-50%)',
                  zIndex: 2,
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
                if (calcItemType(newVal) === ItemTypesEnum.object) {
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
      {isTree && (<>
        {Object.entries(itemValue).sort().map(([k, v]) => {
          const isDir = parentKey === undefined && itemKey === DESC_FIELD && k === FILES_DIRECTORY_FIELD;

          const changeKey = (newKey: string | undefined) => {
            const item = { ...itemValue };
            delete item[k];
            if (newKey !== undefined) {
              item[newKey] = v;
            }
            setItemValue(item);
          };
          return (
            <DetailsViewItem
              key={k}
              itemKey={k}
              itemValue={v}
              parentKey={key}
              setItemValue={newVal => setItemValue({ ...itemValue, [k]: newVal })}
              setItemKey={isDir ? undefined : changeKey}
              def={def?.children?.[k]}
              renderBeforeName={renderBeforeName}
              onRowContextMenu={onRowContextMenu}
              onValueCaret={onValueCaret}
            />
          )
        })}
      </>)}
    </TreeItem>
  )
}
