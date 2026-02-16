import { Stack, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { Toolkit } from "@shared/types";
import { ProjectObj, DocumentObj } from "../../objects/ProjectObj";
import { useProjectStore } from "../../stores/useProjectStore";
import { AddDocumentButton } from "./AddDocumentButton";
import { ProjectDocumentItem } from "./ProjectDocumentItem";
import { compareJsons, filterAndSortByGroups, getValueAtPath } from "../../utils/compareJsons";

const DocumentSplitGroup = ({
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
    const val = getValueAtPath(doc.data.desc, bestPath);
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

export const ToolkitTreeItem = ({
  project,
  toolkit,
  showEmpty = false,
  showDocumentPreview = true,
  maxDepth = 3,
  k = 2,
}: {
  toolkit: Toolkit | undefined;
  project: ProjectObj;
  showEmpty?: boolean;
  showDocumentPreview?: boolean;
  maxDepth?: number;
  k?: number;
}) => {
  const { toolkits } = useProjectStore();

  const toolkitName = toolkit ? toolkit.toolkit : "no-toolkit";
  const toolkitLabel = toolkit ? toolkit.toolkit : "No Toolkit Documents";
  const docs =
    project?.documents.filter((d) => {
      if (toolkit) {
        return d.toolkit === toolkitName;
      }
      return !toolkits.some((t) => t.toolkit === d.toolkit);
    }) || [];

  return !docs.length && !showEmpty ? null : (
    <TreeItem
      key={toolkitName}
      itemId={toolkitName}
      label={
        <Stack direction="row" spacing={1} justifyContent="start" alignItems="center">
          <Typography>{toolkitLabel}</Typography>
          <AddDocumentButton toolkit={toolkit} />
        </Stack>
      }
    >
      <DocumentSplitGroup
        docs={docs}
        project={project}
        showDocumentPreview={showDocumentPreview}
        depth={maxDepth}
        k={k}
      />
    </TreeItem>
  );
};
