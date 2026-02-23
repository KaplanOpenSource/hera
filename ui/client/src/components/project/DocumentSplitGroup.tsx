import { Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { DocumentObj, ProjectObj } from "../../objects/ProjectObj";
import { compareJsons, CompareResult, filterAndSortByGroups, getValueAtPath } from "../../utils/compareJsons";
import { ProjectDocumentItem } from "./ProjectDocumentItem";
import { Case, SwitchCase } from "../../elements/SwitchCase";

const VALUE_GROUP_REST = "__rest__";
const VALUE_GROUP_UNDEFINED = "__undefined__";

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
  minGroupSize,
  compared,
}: {
  docs: DocumentObj[];
  project: ProjectObj;
  showDocumentPreview: boolean;
  depth: number;
  minGroupSize: number;
  compared: CompareResult[];
}) => {
  const bestPath = compared[0].path;
  const fieldLabel = bestPath.replace(/^\//, "").replace(/\//g, ".");

  const valueCounts = new Map<string, number>();
  for (const doc of docs) {
    const val = getValueAtPath(doc.data.desc as any, bestPath);
    const key = val === undefined ? VALUE_GROUP_UNDEFINED : String(val);
    valueCounts.set(key, (valueCounts.get(key) || 0) + 1);
  }

  const groups = new Map<string, DocumentObj[]>();
  const restDocs: DocumentObj[] = [];
  for (const doc of docs) {
    const val = getValueAtPath(doc.data.desc as any, bestPath);
    const key = val === undefined ? VALUE_GROUP_UNDEFINED : String(val);
    if (valueCounts.get(key)! < minGroupSize) {
      restDocs.push(doc);
    } else {
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key)!.push(doc);
    }
  }

  if (restDocs.length > 0) {
    groups.set(VALUE_GROUP_REST, restDocs);
  }

  return (
    <>
      {[...groups.entries()].map(([value, groupDocs]) => {
        const itemKey = `split_${fieldLabel}=${value}`;
        return (
          <TreeItem
            key={itemKey}
            itemId={itemKey}
            label={
              <Typography>
                <SwitchCase test={value}>
                  <Case value={VALUE_GROUP_REST}>
                    <b>{fieldLabel}</b> other values
                  </Case>
                  <Case value={VALUE_GROUP_UNDEFINED}>
                    <b>{fieldLabel}</b> not existing
                  </Case>
                  <Case isDefault={true}>
                    <b>{fieldLabel}</b> == <b>{value}</b>
                  </Case>
                </SwitchCase>
              </Typography>
            }
          >
            <DocumentSplitGroup
              docs={groupDocs}
              project={project}
              showDocumentPreview={showDocumentPreview}
              minGroupSize={minGroupSize}
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
  minGroupSize,
}: {
  docs: DocumentObj[];
  project: ProjectObj;
  showDocumentPreview: boolean;
  depth: number;
  minGroupSize: number;
}) => {
  let compared: CompareResult[] = [];
  if (depth > 0 && docs.length > 1) {
    const descs = docs.map((d) => d.data.desc);
    compared = filterAndSortByGroups(compareJsons(descs, true), minGroupSize);
  }

  return compared.length
    ? (
      <DocumentSplitTree
        docs={docs}
        project={project}
        showDocumentPreview={showDocumentPreview}
        depth={depth}
        minGroupSize={minGroupSize}
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