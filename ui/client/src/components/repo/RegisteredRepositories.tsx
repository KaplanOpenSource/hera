import { TreeItem } from "@mui/x-tree-view";
import { RepoContents } from "./RepoContents";
import { useRegisteredRepositories } from "./useRegisteredRepositories";
import { DEFAULT_PROJECT } from "../../stores/useProjectStore";

export const RegisteredRepositories = () => {
  const { repositories } = useRegisteredRepositories();

  return (
    <TreeItem itemId="registered-repos" label={<>
      Repositories in <b>{DEFAULT_PROJECT}</b>
    </>}>
      {repositories.length === 0
        ? <TreeItem itemId="registered-repos/__empty" label="No registered repositories" />
        : repositories.map(r => (
          <RepoContents key={r.datasourceName} repo={r} />
        ))
      }
    </TreeItem>
  )
}
