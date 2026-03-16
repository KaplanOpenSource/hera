import { DocumentObj } from "../objects/ProjectObj";
import { ViewSettingsType } from "../stores/useViewSettingsStore";
import { compareJsons, CompareResult, filterAndSortByGroups, getValueAtPath } from "./compareJsons";
import { reorderByKeys } from "./utils";

export const VALUE_GROUP_REST = "__rest__";
export const VALUE_GROUP_UNDEFINED = "__undefined__";
export const DESC_PATH_TOOLKIT = "/toolkit";
export const DESC_PATH_TYPE = "/type";

export enum SplitTreeNodeType {
  Split = 'split',
  Leaf = 'leaf',
}

export type SplitNode = {
  type: SplitTreeNodeType.Split;
  itemKey: string;
  path: string;
  value: string;
  children: SplitTreeNode[];
};

export type LeafNode = {
  type: SplitTreeNodeType.Leaf;
  doc: DocumentObj;
};

export type SplitTreeNode = SplitNode | LeafNode;

function buildSplitLevel(
  docs: DocumentObj[],
  compared: CompareResult[],
  viewSettings: ViewSettingsType,
): { itemKey: string; value: string; path: string; docs: DocumentObj[] }[] {
  const bestPath = compared[0].path;
  const isToolkit = bestPath === DESC_PATH_TOOLKIT;

  const valueCounts = new Map<string, number>();
  for (const doc of docs) {
    const val = getValueAtPath(doc.data.desc as any, bestPath);
    const key = val === undefined ? VALUE_GROUP_UNDEFINED : String(val);
    valueCounts.set(key, (valueCounts.get(key) || 0) + 1);
  }

  const groups = new Map<string, DocumentObj[]>();
  const restDocs: DocumentObj[] = [];
  for (const doc of docs) {
    const val = getValueAtPath(doc.extDesc, bestPath);
    const key = val === undefined ? VALUE_GROUP_UNDEFINED : String(val);
    if (!isToolkit && valueCounts.get(key)! < viewSettings.minGroupSize) {
      restDocs.push(doc);
    } else {
      if (!groups.has(key)) {
        groups.set(key, []);
      }
      groups.get(key)!.push(doc);
    }
  }

  if (restDocs.length > 0) {
    groups.set(VALUE_GROUP_REST, restDocs);
  }

  return [...groups.entries()]
    .sort((a, b) => a[0].localeCompare(b[0]))
    .map(([value, groupDocs]) => ({
      itemKey: `split_${bestPath}=${value}`,
      value,
      path: bestPath,
      docs: groupDocs,
    }));
}

function getCompared(docs: DocumentObj[], depth: number, viewSettings: ViewSettingsType): CompareResult[] {
  if (depth <= 0 || docs.length <= 1) return [];
  const descs = docs.map((d) => ({ ...d.extDesc }));
  const paths = compareJsons(descs, true);
  let compared = filterAndSortByGroups(paths, viewSettings.minGroupSize, viewSettings.maxBranches);
  if (viewSettings.firstBranchHeadFields) {
    compared = reorderByKeys(compared, x => x.path, [DESC_PATH_TOOLKIT, DESC_PATH_TYPE]);
  }
  return compared;
}

export function buildSplitTree(
  docs: DocumentObj[],
  depth: number,
  viewSettings: ViewSettingsType,
): SplitTreeNode[] {
  const compared = getCompared(docs, depth, viewSettings);
  if (!compared.length) {
    return docs.map(doc => ({ type: SplitTreeNodeType.Leaf, doc }));
  }
  const groups = buildSplitLevel(docs, compared, viewSettings);
  return groups.map(g => ({
    type: SplitTreeNodeType.Split,
    itemKey: g.itemKey,
    path: g.path,
    value: g.value,
    children: buildSplitTree(g.docs, depth - 1, viewSettings),
  }));
}
