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
  return (
    <TreeItem
      key={key}
      itemId={key}
      label={itemKey}
    >
      {typeof itemValue === 'object'
        ? Object.entries(itemValue).map(([k, v], i) => (
          <DetailsViewItem
            key={k}
            itemKey={k}
            itemValue={v}
            level={level + 1}
            index={i}
          />
        ))
        : (
          <span>
            {JSON.stringify(itemValue)}
          </span>
        )}
    </TreeItem>

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