import { SimpleTreeView } from "@mui/x-tree-view";
import { RepoContents } from "./RepoContents";
import { useRegisteredRepositories } from "./useRegisteredRepositories";

export const RepositoriesInProject = () => {
  const { repositories } = useRegisteredRepositories();

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
