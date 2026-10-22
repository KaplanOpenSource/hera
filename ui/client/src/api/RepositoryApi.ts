import { Repository } from "../shared/types";
import { ApiCall } from "./ApiCall";

export type RepositoryTrees = {
  // The repository file as it is on disk.
  json: any,
  // The same tree, with the datasource paths made absolute.
  resolved: any,
};

export type RepoJsonsInFolder = {
  // The repository files found under the folder, sorted.
  files: string[],
  // Path to contents, empty unless the contents were asked for.
  repoJsons: { [path: string]: any },
};

// Every call to /repository/*. Each method answers data, or null when the call
// failed.
export class RepositoryApi extends ApiCall {
  // One repository file, by its path. A repository name never crosses the wire.
  async get(path: string): Promise<RepositoryTrees | null> {
    const data = await this.getJson('/repository/get', { path }, { label: 'repo JSON' });
    if (!data) {
      return null;
    }
    return { json: data.json, resolved: data.resolved };
  }

  // The repositories registered in the default project, for showing their names.
  async registered(): Promise<Repository[]> {
    const data = await this.getJson('/repository/registered', {}, { label: 'registered repositories' });
    if (!data) {
      return [];
    }
    return (data.repositories ?? []) as Repository[];
  }

  // The repository files under a folder. A folder that is not there is empty.
  async jsonsInFolder(folder: string, withContents = false): Promise<RepoJsonsInFolder | null> {
    const data = await this.getJson('/repository/jsonsInFolder',
      { folder, withContents: String(withContents) },
      { label: 'central repo files' });
    if (!data) {
      return null;
    }
    return { files: (data.files ?? []) as string[], repoJsons: data.repoJsons ?? {} };
  }

  // One call for all three sources. Its answer carries the project, so nothing reloads it.
  async loadIntoProject(body: object, label: string): Promise<any | null> {
    return await this.postJson('/repository/loadIntoProject', body, { label });
  }
}
