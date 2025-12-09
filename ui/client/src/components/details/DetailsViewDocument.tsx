import { Stack, TextField, Typography } from '@mui/material';
import { SimpleTreeView, TreeItem } from '@mui/x-tree-view';

const DetailsViewItem = ({
  itemKey,
  itemValue,
  level = 0,
  index,
}: {
  itemKey: string,
  itemValue: any,
  level: number,
  index: number,
}) => {
  const key = `___lvl${level}_idx${index}_${itemKey}`
  return typeof itemValue === 'object'
    ? (
      <TreeItem
        key={key}
        itemId={key}
        label={itemKey}
      >
        {Object.entries(itemValue).map(([k, v], i) => (
          <DetailsViewItem
            key={k}
            itemKey={k}
            itemValue={v}
            level={level + 1}
            index={i}
          />
        ))}
      </TreeItem>
    )
    : (
      <TreeItem
        key={key}
        itemId={key}
        label={(
          <Stack direction='row' spacing={1} justifyItems={'center'} alignItems={'center'}>
            <Typography>
              {itemKey}
            </Typography>
            <TextField
              size='small'
              value={
                JSON.stringify(itemValue)
              }
            />
          </Stack>
        )}
      />
    )
}

export const DetailsViewDocument = ({
  doc,
}: {
  doc: any,
}) => {
  const name = doc?.desc?.datasourceName || doc?.type || doc._cls;
  const items: [string, any][] = Object.entries(doc);

  return (
    <SimpleTreeView>
      <DetailsViewItem
        itemKey={name}
        itemValue={doc}
        level={0}
        index={0}
      />
    </SimpleTreeView>
  )
}