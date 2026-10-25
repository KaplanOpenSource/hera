import { Typography } from "@mui/material";
import { useEffect, useState } from "react";
import { repositoryApi } from "../../api/api";
import { RepoJsonMerger } from "./RepoJsonMerger";
import { RepoTreeDisplay } from "./RepoTreeDisplay";

const STORAGE_KEY = 'hera-central-repo-folder';
const DEFAULT_FOLDER = '~/hera/repositories/';

export const DetailsViewMergedRepo = () => {
  const [merged, setMerged] = useState<{ [key: string]: any } | null>(null);
  const [overrides, setOverrides] = useState<{ [path: string]: string[] }>({});
  const [repoJsons, setRepoJsons] = useState<{ [path: string]: { [key: string]: any } }>({});
  const [loading, setLoading] = useState(true);

  const folder = localStorage.getItem(STORAGE_KEY) || DEFAULT_FOLDER;

  useEffect(() => {
    (async () => {
      setLoading(true);
      const answer = await repositoryApi.jsonsInFolder(folder, true);
      if (answer) {
        const merger = new RepoJsonMerger(answer.repoJsons);
        setMerged(merger.merged);
        setOverrides(merger.overrides);
        setRepoJsons(answer.repoJsons);
      }
      setLoading(false);
    })();
  }, [folder]);

  if (loading) {
    return <Typography>Loading...</Typography>;
  }

  if (!merged || Object.keys(merged).length === 0) {
    return <Typography>No valid repository files found in {folder}</Typography>;
  }

  return (
    <RepoTreeDisplay
      tree={merged}
      setTree={setMerged}
      label={`${folder} (merged)`}
      itemId="merged-root"
      overrides={overrides}
      repoJsons={repoJsons}
    />
  );
};
