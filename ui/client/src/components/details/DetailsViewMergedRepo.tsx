import { Typography } from "@mui/material";
import { useEffect, useState } from "react";
import { fetchPython } from "../../io/fetchPython";
import { LOAD_REPO_JSON_PYTHON } from "../../shared/repoJsonPython";
import { RepoTreeDisplay } from "./RepoTreeDisplay";

const STORAGE_KEY = 'hera-central-repo-folder';
const DEFAULT_FOLDER = '~/hera/repositories/';

const mergeRepoDocs = (docs: Record<string, any>[]) => {
  const merged: Record<string, any> = {};
  for (const doc of docs) {
    for (const [toolkit, sections] of Object.entries(doc)) {
      if (!merged[toolkit]) {
        merged[toolkit] = {};
      }
      if (toolkit.toLowerCase() === 'experiment') {
        const dsSection = Object.entries(sections as Record<string, any>)
          .find(([key]) => key.toLowerCase() === 'datasource');
        const dsName = dsSection
          && typeof dsSection[1] === 'object'
          && dsSection[1] !== null
          && Object.keys(dsSection[1]).length > 0
          ? Object.keys(dsSection[1])[0]
          : null;
        if (dsName) {
          if (!merged[toolkit][dsName]) {
            merged[toolkit][dsName] = {};
          }
          for (const [section, sectionData] of Object.entries(sections as Record<string, any>)) {
            merged[toolkit][dsName][section] = {
              ...merged[toolkit][dsName][section],
              ...(sectionData as Record<string, any>),
            };
          }
          continue;
        }
      }
      for (const [section, sectionData] of Object.entries(sections as Record<string, any>)) {
        merged[toolkit][section] = {
          ...merged[toolkit][section],
          ...(sectionData as Record<string, any>),
        };
      }
    }
  }
  return merged;
};

export const DetailsViewMergedRepo = () => {
  const [merged, setMerged] = useState<Record<string, any> | null>(null);
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
repoJsons = []
if os.path.isdir(folder):
    allFiles = glob.glob(os.path.join(folder, '**', '*.json'), recursive=True)
    allFiles = [f for f in allFiles if not f.endswith('caseConfiguration.json')]
    for f in sorted(allFiles):
        doc = loadRepoJson(f)
        if doc is not None:
            repoJsons.append(doc)
`,
      });
      if (data?.repoJsons) {
        setMerged(mergeRepoDocs(data.repoJsons));
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
