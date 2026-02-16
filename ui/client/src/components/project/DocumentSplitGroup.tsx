import { Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { DocumentObj, ProjectObj } from "../../objects/ProjectObj";
import { compareJsons, filterAndSortByGroups, getValueAtPath } from "../../utils/compareJsons";
import { ProjectDocumentItem } from "./ProjectDocumentItem";

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
  if (depth <= 0 || docs.length <= 1) {
    return (
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
  }

  const compared = compareJsons(docs.map((d) => d.data.desc), true);
  const sorted = filterAndSortByGroups(compared, k);

  if (!sorted.length) {
    return (
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
  }

  const bestPath = sorted[0].path;

  const groups = new Map<string, DocumentObj[]>();
  for (const doc of docs) {
    const val = getValueAtPath(doc.data.desc as any, bestPath);
    const key = val === undefined ? "__undefined__" : String(val);
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key)!.push(doc);
  }

  return (
    <>
      {[...groups.entries()].map(([value, groupDocs]) => (
        <TreeItem
          key={`${bestPath}=${value}`}
          itemId={`split_${bestPath}=${value}`}
          label={
            <Typography variant="body2">
              {bestPath}: {value === "__undefined__" ? "—" : value}
            </Typography>
          }
        >
          <DocumentSplitGroup
            docs={groupDocs}
            project={project}
            showDocumentPreview={showDocumentPreview}
            depth={depth - 1}
            k={k}
          />
        </TreeItem>
      ))}
    </>
  );
};
