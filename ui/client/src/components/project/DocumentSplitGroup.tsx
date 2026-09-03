import { TreeItem } from "@mui/x-tree-view";
import { DocumentObj, ProjectObj } from "../../objects/ProjectObj";
import { useViewSettingsStore } from "../../stores/useViewSettingsStore";
import { classifyDocument } from "../../shared/tabKind";
import { SplitTree, SplitTreeNode, SplitTreeNodeType } from "../../utils/splitTree";
import { DocumentSplitTreeLabel } from "./DocumentSplitTreeLabel";
import { ProjectDocumentItem } from "./ProjectDocumentItem";

const DocumentSplitTree = ({
  project,
  nodes,
  onDocumentDeleted,
}: {
  project: ProjectObj;
  nodes: SplitTreeNode[];
  onDocumentDeleted?: () => void;
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
              kind={classifyDocument(node.doc)}
              onDocumentDeleted={onDocumentDeleted}
            />
          );
        }
        // console.log('[tree-item] split:', node.itemKey);
        return (
          <TreeItem
            key={node.itemKey}
            itemId={node.itemKey}
            label={
              <DocumentSplitTreeLabel
                itemKey={node.itemKey}
                path={node.path}
                value={node.value}
              />
            }
          >
            <DocumentSplitTree
              project={project}
              nodes={node.children}
              onDocumentDeleted={onDocumentDeleted}
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
  onDocumentDeleted,
}: {
  docs: DocumentObj[];
  project: ProjectObj;
  depth: number;
  onDocumentDeleted?: () => void;
}) => {
  const { viewSettings } = useViewSettingsStore();
  const tree = new SplitTree(docs, depth, viewSettings);
  return (
    <DocumentSplitTree
      project={project}
      nodes={tree.nodes}
      onDocumentDeleted={onDocumentDeleted}
    />
  );
};
