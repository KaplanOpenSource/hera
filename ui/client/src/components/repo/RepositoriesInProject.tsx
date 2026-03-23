import { SimpleTreeView } from "@mui/x-tree-view";
import { Repository } from "../../shared/types";
import { RepoContents } from "./RepoContents";

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
        <RepoContents key={r.datasourceName} repo={r} />
      ))}
    </SimpleTreeView>
  )
}
