import { Link, LinkOff } from "@mui/icons-material";
import { useEffect, useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { repositoryApi } from "../../api/api";
import { idRepoId, TEMP_REPO_NAME } from "../../shared/idDocId";
import { RepoTreeDisplay } from "./RepoTreeDisplay";

export const DetailsViewRepo = ({
  repoPath,
}: {
  repoPath: string,
}) => {
  const isTempRepo = repoPath.includes(TEMP_REPO_NAME);
  const [tree, setTree] = useState<any>(undefined);
  const [origTree, setOrigTree] = useState<any>(undefined);
  const [peekOriginal, setPeekOriginal] = useState(false);

  useEffect(() => {
    (async () => {
      if (!isTempRepo) {
        const trees = await repositoryApi.get(repoPath);
        setOrigTree(trees?.json ?? undefined);
        setTree(trees?.resolved ?? trees?.json ?? undefined);
      }
    })()
  }, [repoPath])

  const pathsResolved = JSON.stringify(origTree) !== JSON.stringify(tree);
  const displayedTree = (peekOriginal && origTree) ? origTree : tree;

  return (
    <RepoTreeDisplay
      tree={displayedTree}
      setTree={setTree}
      label={repoPath.split('/').pop() ?? repoPath}
      itemId={idRepoId(repoPath)}
      showStrDefault={isTempRepo}
      labelIcons={pathsResolved && (
        <ButtonTooltip
          title={peekOriginal ? 'Original (relative) paths' : 'Resolved (absolute) paths — hover to see original'}
          onClick={() => { }}
          onMouseEnter={() => setPeekOriginal(true)}
          onMouseLeave={() => setPeekOriginal(false)}
        >
          {peekOriginal ? <LinkOff fontSize="small" /> : <Link fontSize="small" />}
        </ButtonTooltip>
      )}
    />
  );
};
