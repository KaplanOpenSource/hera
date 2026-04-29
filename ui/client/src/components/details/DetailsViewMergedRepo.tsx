import { Typography } from "@mui/material";
import { useEffect, useState } from "react";
import { fetchPython } from "../../io/fetchPython";
import { LOAD_REPO_JSON_PYTHON } from "../../shared/repoJsonPython";
import { RepoTreeDisplay } from "./RepoTreeDisplay";

const STORAGE_KEY = 'hera-central-repo-folder';
const DEFAULT_FOLDER = '~/hera/repositories/';

export const DetailsViewMergedRepo = () => {
  const [merged, setMerged] = useState<Record<string, any> | null>(null);
  const [loading, setLoading] = useState(true);

  const folder = localStorage.getItem(STORAGE_KEY) || DEFAULT_FOLDER;

  useEffect(() => {
    (async () => {
      setLoading(true);
      const { data } = await fetchPython({
        results: ['mergedJson'],
        label: 'merged central repo',
        code: `
import os, glob, json
${LOAD_REPO_JSON_PYTHON}
folder = os.path.expanduser('${folder}')
merged = {}
if os.path.isdir(folder):
    allFiles = glob.glob(os.path.join(folder, '**', '*.json'), recursive=True)
    allFiles = [f for f in allFiles if not f.endswith('caseConfiguration.json')]
    for f in sorted(allFiles):
        doc = loadRepoJson(f)
        if doc is not None:
            for toolkit, sections in doc.items():
                if toolkit not in merged:
                    merged[toolkit] = {}
                for section, sectionData in sections.items():
                    if section not in merged[toolkit]:
                        merged[toolkit][section] = {}
                    merged[toolkit][section].update(sectionData)
mergedJson = merged
`,
      });
      if (data?.mergedJson) {
        setMerged(data.mergedJson);
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
    />
  );
};
