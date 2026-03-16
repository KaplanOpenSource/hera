import { DocumentObj } from "../objects/ProjectObj";
import { ViewSettingsType } from "../stores/useViewSettingsStore";
import { compareJsons, filterAndSortByGroups, getValueAtPath } from "./compareJsons";
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

export class SplitTree {
  public readonly nodes: SplitTreeNode[];
  private readonly docIds: Set<string>;

  constructor(
    public readonly docs: DocumentObj[],
    public readonly depth: number,
    public readonly viewSettings: ViewSettingsType,
  ) {
    this.nodes = this.buildNodes(docs, depth);
    this.docIds = new Set(docs.map(d => d.docid));
  }

  findAncestorKeys(docId: string): string[] | null {
    return this.findAncestors(this.nodes, docId);
  }

  needsRebuild(docs: DocumentObj[], viewSettings: ViewSettingsType): boolean {
    if (viewSettings !== this.viewSettings) return true;
    if (docs.length !== this.docIds.size) return true;
    return docs.some(d => !this.docIds.has(d.docid));
  }

  private buildNodes(docs: DocumentObj[], depth: number): SplitTreeNode[] {
    const compared = this.getCompared(docs, depth);
    if (!compared.length) {
      return docs.map(doc => ({ type: SplitTreeNodeType.Leaf, doc }));
    }
    const groups = this.buildSplitLevel(docs, compared);
    return groups.map(g => ({
      type: SplitTreeNodeType.Split,
      itemKey: g.itemKey,
      path: g.path,
      value: g.value,
      children: this.buildNodes(g.docs, depth - 1),
    }));
  }

  private getCompared(docs: DocumentObj[], depth: number) {
    if (depth <= 0 || docs.length <= 1) return [];
    const descs = docs.map((d) => ({ ...d.extDesc }));
    const paths = compareJsons(descs, true);
    let compared = filterAndSortByGroups(paths, this.viewSettings.minGroupSize, this.viewSettings.maxBranches);
    if (this.viewSettings.firstBranchHeadFields) {
      compared = reorderByKeys(compared, x => x.path, [DESC_PATH_TOOLKIT, DESC_PATH_TYPE]);
    }
    return compared;
  }

  private buildSplitLevel(docs: DocumentObj[], compared: { path: string }[]) {
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
      if (!isToolkit && valueCounts.get(key)! < this.viewSettings.minGroupSize) {
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

  private findAncestors(nodes: SplitTreeNode[], docId: string): string[] | null {
    for (const node of nodes) {
      if (node.type === SplitTreeNodeType.Leaf) {
        if (node.doc.docid === docId) return [];
      } else {
        const result = this.findAncestors(node.children, docId);
        if (result !== null) return [node.itemKey, ...result];
      }
    }
    return null;
  }
}
