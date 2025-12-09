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
      ret.push({ id: k, label: k, children: [{ id: k + '_ch', label: JSON.stringify(v) }] });
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