import { SimpleTreeView, TreeItem } from "@mui/x-tree-view";
import { Repository } from "../../shared/types";

export const RepositoriesInProject = ({
  repositories,
}: {
  repositories: Repository[],
}) => {
  if (repositories.length === 0) {
    return null;
  }
  return (
    <SimpleTreeView>
      {repositories.map(r => (
        <TreeItem key={r.datasourceName} itemId={r.datasourceName} label={r.datasourceName}>
          <TreeItem itemId={`${r.datasourceName}/resource`} label={`Resource: ${r.resource}`} />
          <TreeItem itemId={`${r.datasourceName}/dataFormat`} label={`Data Format: ${r.dataFormat}`} />
          <TreeItem itemId={`${r.datasourceName}/toolkit`} label={`Toolkit: ${r.toolkit}`} />
          <TreeItem itemId={`${r.datasourceName}/version`} label={`Version: ${r.version.join('.')}`} />
        </TreeItem>
      ))}
    </SimpleTreeView>
  )
}
