import { TreeItem } from "@mui/x-tree-view";
import { DocumentObj, ProjectObj } from "../../objects/ProjectObj";
import { useViewSettingsStore } from "../../stores/useViewSettingsStore";
import { compareJsons, CompareResult, filterAndSortByGroups, getValueAtPath } from "../../utils/compareJsons";
import { DocumentSplitTreeLabel } from "./DocumentSplitTreeLabel";
import { ProjectDocumentItem } from "./ProjectDocumentItem";

export const VALUE_GROUP_REST = "__rest__";
export const VALUE_GROUP_UNDEFINED = "__undefined__";
export const TOOLKIT_DESC_PATH = "/toolkit";

const DocumentSplitTree = ({
  docs,
  project,
  depth,
  compared,
}: {
  docs: DocumentObj[];
  project: ProjectObj;
  depth: number;
  compared: CompareResult[];
}) => {
  const { viewSettings } = useViewSettingsStore();
  const bestPath = compared[0].path;
  const isToolkit = bestPath === TOOLKIT_DESC_PATH;

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

  const entries = [...groups.entries()].sort((a, b) => (a[0].localeCompare(b[0])));

  return (
    <>
      {entries.map(([value, groupDocs]) => {
        const itemKey = `split_${bestPath}=${value}`;
        return (
          <TreeItem
            key={itemKey}
            itemId={itemKey}
            label={
              <DocumentSplitTreeLabel
                path={bestPath}
                value={value}
              />
            }
          >
            <DocumentSplitGroup
              docs={groupDocs}
              project={project}
              depth={depth - 1}
            />
          </TreeItem>
        );
      })}
    </>
  );
};

export const DocumentSplitGroup = ({
  docs,
  project,
  depth,
}: {
  docs: DocumentObj[];
  project: ProjectObj;
  depth: number;
}) => {
  const { viewSettings } = useViewSettingsStore();
  let compared: CompareResult[] = [];
  if (depth > 0 && docs.length > 1) {
    const descs = docs.map((d) => ({ ...d.extDesc }));
    const paths = compareJsons(descs, true);
    compared = filterAndSortByGroups(paths, viewSettings.minGroupSize, viewSettings.maxBranches);
    if (viewSettings.firstBranchByToolkits) {
      const i = compared.findIndex(x => x.path === TOOLKIT_DESC_PATH);
      if (i > 0) {
        compared = [compared[i], ...compared.filter((_, j) => i !== j)];
      }
    }
  }

  return compared.length
    ? (
      <DocumentSplitTree
        docs={docs}
        project={project}
        depth={depth}
        compared={compared}
      />
    )
    : (
      <>
        {docs.map((doc) => (
          <ProjectDocumentItem
            key={`proj${project.name}_doc${doc.docid}`}
            project={project.data}
            document={doc.data}
          />
        ))}
      </>
    );
};