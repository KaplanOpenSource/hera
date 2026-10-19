import { getJson } from "./entrypoint";

export type RepositoryTrees = {
  // The repository file as it is on disk.
  json: any,
  // The same tree, with the datasource paths made absolute.
  resolved: any,
};

// One repository file, by its path. A repository name never crosses the wire.
export const fetchRepository = async (path: string): Promise<RepositoryTrees | null> => {
  const data = await getJson('/repository/get', { path }, { label: 'repo JSON' });
  if (!data) {
    return null;
  }
  return { json: data.json, resolved: data.resolved };
};
