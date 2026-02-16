import { Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { DocumentObj, ProjectObj } from "../../objects/ProjectObj";
import { compareJsons, CompareResult, filterAndSortByGroups, getValueAtPath } from "../../utils/compareJsons";
import { ProjectDocumentItem } from "./ProjectDocumentItem";

const DocumentItems = ({
  docs,
  project,
  showDocumentPreview,
}: {
  docs: DocumentObj[];
  project: ProjectObj;
  showDocumentPreview: boolean;
}) => (
  <>
    {docs.map((doc) => (
      <ProjectDocumentItem
        key={`proj${project.name}_doc${doc.docid}`}
        project={project.data}
        document={doc.data}
        showDocumentPreview={showDocumentPreview}
      />
    ))}
  </>
);

const DocumentSplitTree = ({
  docs,
  project,
  showDocumentPreview,
  depth,
  k,
  compared,
}: {
  docs: DocumentObj[];
  project: ProjectObj;
  showDocumentPreview: boolean;
  depth: number;
  k: number;
  compared: CompareResult[];
}) => {
  const bestPath = compared[0].path;
  const fieldLabel = bestPath.replace(/^\//, "").replace(/\//g, ".");

  const valueCounts = new Map<string, number>();
  for (const doc of docs) {
    const val = getValueAtPath(doc.data.desc as any, bestPath);
    const key = val === undefined ? "__undefined__" : String(val);
    valueCounts.set(key, (valueCounts.get(key) || 0) + 1);
  }

  const groups = new Map<string, DocumentObj[]>();
  const restDocs: DocumentObj[] = [];
  for (const doc of docs) {
    const val = getValueAtPath(doc.data.desc as any, bestPath);
    const key = val === undefined ? "__undefined__" : String(val);
    if (valueCounts.get(key)! < k) {
      restDocs.push(doc);
    } else {
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key)!.push(doc);
    }
  }

  if (restDocs.length > 0) {
    groups.set("__rest__", restDocs);
  }

  return (
    <>
      {[...groups.entries()].map(([value, groupDocs]) => {
        const isRest = value === "__rest__";
        const displayValue = isRest
          ? `Documents with other ${fieldLabel} values`
          : `Documents with ${fieldLabel} == ${value === "__undefined__" ? "—" : value}`;
        const itemKey = `split_${fieldLabel}=${value}`;
        return (
          <TreeItem
            key={itemKey}
            itemId={itemKey}
            label={
              <Typography>
                {displayValue}
              </Typography>
            }
          >
            <DocumentSplitGroup
              docs={groupDocs}
              project={project}
              showDocumentPreview={showDocumentPreview}
              k={k}
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
  showDocumentPreview,
  depth,
  k,
}: {
  docs: DocumentObj[];
  project: ProjectObj;
  showDocumentPreview: boolean;
  depth: number;
  k: number;
}) => {
  let compared: CompareResult[] = [];
  if (depth > 0 && docs.length > 1) {
    const descs = docs.map((d) => d.data.desc);
    compared = filterAndSortByGroups(compareJsons(descs, true), k);
  }

  return compared.length
    ? (
      <DocumentSplitTree
        docs={docs}
        project={project}
        showDocumentPreview={showDocumentPreview}
        depth={depth}
        k={k}
        compared={compared}
      />
    )
    : (
      <DocumentItems
        docs={docs}
        project={project}
        showDocumentPreview={showDocumentPreview}
      />
    );
};