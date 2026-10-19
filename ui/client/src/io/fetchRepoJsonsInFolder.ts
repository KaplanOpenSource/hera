import { getJson } from "./entrypoint";

export type RepoJsonsInFolder = {
  // The repository files found under the folder, sorted.
  files: string[],
  // Path to contents, empty unless the contents were asked for.
  repoJsons: { [path: string]: any },
};

// The repository files under a folder. A folder that is not there is empty.
export const fetchRepoJsonsInFolder = async (
  folder: string,
  withContents = false,
): Promise<RepoJsonsInFolder | null> => {
  const data = await getJson('/repository/jsonsInFolder',
    { folder, withContents: String(withContents) },
    { label: 'central repo files' });
  if (!data) {
    return null;
  }
  return { files: (data.files ?? []) as string[], repoJsons: data.repoJsons ?? {} };
};
