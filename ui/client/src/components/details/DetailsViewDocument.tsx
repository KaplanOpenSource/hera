import { RichTreeView } from '@mui/x-tree-view';
import { TreeViewBaseItem } from '@mui/x-tree-view/models';

export const DetailsViewDocument = ({
  doc,
}: {
  doc: any,
}) => {
  const name = doc?.desc?.datasourceName || doc?.type || doc._cls;

  const convert = (s: any) => {
    const ret: TreeViewBaseItem[] = [];
    for (const [k, v] of Object.entries(s)) {
      const children = typeof v === 'object' ? convert(v) : [{ id: k + '_ch', label: JSON.stringify(v) }];
      ret.push({ id: k, label: k, children });
    }
    return ret;
  }
  const obj: TreeViewBaseItem = {
    id: '____rooot____',
    label: name,
    children: convert(doc)
  };

  return (<RichTreeView
    items={[obj]}
  />)
}