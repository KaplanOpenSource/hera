import { Typography } from "@mui/material";
import { useEffect, useState } from "react";
import { fetchPython } from "../../io/fetchPython";
import { LOAD_REPO_JSON_PYTHON } from "../../shared/repoJsonPython";
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
      const { data } = await fetchPython({
        results: ['repoJsons'],
        label: 'merged central repo',
        code: `
import os, glob, json
${LOAD_REPO_JSON_PYTHON}
folder = os.path.expanduser('${folder}')
repoJsons = {}
if os.path.isdir(folder):
    allFiles = glob.glob(os.path.join(folder, '**', '*.json'), recursive=True)
    allFiles = [f for f in allFiles if not f.endswith('caseConfiguration.json')]
    for f in sorted(allFiles):
        doc = loadRepoJson(f)
        if doc is not None:
            repoJsons[f] = doc
`,
      });
      if (data?.repoJsons) {
        const merger = new RepoJsonMerger(data.repoJsons);
        setMerged(merger.merged);
        setOverrides(merger.overrides);
        setRepoJsons(data.repoJsons);
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
