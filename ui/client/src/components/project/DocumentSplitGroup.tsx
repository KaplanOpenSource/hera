import { TreeItem } from "@mui/x-tree-view";
import { DocumentObj, ProjectObj } from "../../objects/ProjectObj";
import { useViewSettingsStore } from "../../stores/useViewSettingsStore";
import { buildSplitTree, SplitTreeNode, SplitTreeNodeType } from "../../utils/splitTree";
import { DocumentSplitTreeLabel } from "./DocumentSplitTreeLabel";
import { ProjectDocumentItem } from "./ProjectDocumentItem";

const DocumentSplitTree = ({
  project,
  nodes,
}: {
  project: ProjectObj;
  nodes: SplitTreeNode[];
}) => {
  return (
    <>
      {nodes.map(node => {
        if (node.type === SplitTreeNodeType.Leaf) {
          return (
            <ProjectDocumentItem
              key={`proj${project.name}_doc${node.doc.docid}`}
              project={project.data}
              document={node.doc.data}
            />
          );
        }
        console.log('[tree-item] split:', node.itemKey);
        return (
          <TreeItem
            key={node.itemKey}
            itemId={node.itemKey}
            label={
              <DocumentSplitTreeLabel
                path={node.path}
                value={node.value}
              />
            }
          >
            <DocumentSplitTree
              project={project}
              nodes={node.children}
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
  const nodes = buildSplitTree(docs, depth, viewSettings);
  return (
    <DocumentSplitTree
      project={project}
      nodes={nodes}
    />
  );
};
