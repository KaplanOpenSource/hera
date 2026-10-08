import { idDocId } from "../shared/idDocId";
import { SplitTreeNode, SplitTreeNodeType } from "./splitTree";

// Finds a branch node by its item key, anywhere in the tree.
export const findNodeByKey = (nodes: SplitTreeNode[], itemKey: string): SplitTreeNode | undefined => {
  for (const node of nodes) {
    if (node.type === SplitTreeNodeType.Split) {
      if (node.itemKey === itemKey) {
        return node;
      }
      const found = findNodeByKey(node.children, itemKey);
      if (found) {
        return found;
      }
    }
  }
  return undefined;
};

// Every tree item key in a subtree: the branches and the documents below them.
export const collectSubtreeKeys = (node: SplitTreeNode): string[] => {
  if (node.type === SplitTreeNodeType.Leaf) {
    return [idDocId(node.doc.docid)];
  }
  const keys = [node.itemKey];
  for (const child of node.children) {
    keys.push(...collectSubtreeKeys(child));
  }
  return keys;
};
