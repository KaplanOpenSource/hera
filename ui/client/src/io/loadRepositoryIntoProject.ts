import { ProjectEntire } from "../shared/types";
import { useProjectStore } from "../stores/useProjectStore";
import { postJson } from "./entrypoint";
import { pushError, pushInfo } from "./snackbar";

type LoadRepositoryBody = {
  projectName: string,
  // A repository file, or the tree itself with the directory its paths are against.
  path?: string,
  repositoryJSON?: object,
  baseDir?: string,
  overwrite: boolean,
};

// hera skips a toolkit it cannot build and a section that fails, so a load can
// add nothing and still succeed. The answer says what happened; this shows it.
const report = (label: string, added: number, problems: string[]) => {
  if (problems.length) {
    console.error(`${label}: ${problems.join('\n')}`);
  }
  if (added === 0) {
    if (problems.length) {
      pushError(`${label}: nothing was added - ${problems[0]}`);
      return;
    }
    // Never guess a reason: say there is none and where the rest would be.
    pushError(`${label}: nothing was added, and hera reported nothing - see the server console`);
    return;
  }
  if (problems.length) {
    pushInfo(`${label}: added ${added}, ${problems.length} skipped - see the console`);
  }
};

// One call for all three sources. Its answer carries the project, so nothing reloads it.
const load = async (body: LoadRepositoryBody, label: string) => {
  const data = await postJson('/repository/loadIntoProject', body, { label });
  if (!data) {
    return false;
  }
  const { currProjectName, setCurrentProject } = useProjectStore.getState();
  // Skip if the user switched projects while this call was in flight.
  if (currProjectName === body.projectName) {
    setCurrentProject(data.project as ProjectEntire);
  }
  report(label, data.added ?? 0, data.problems ?? []);
  return true;
};

// Loads the datasources of one repository file into the current project.
export const loadRepositoryIntoProject = async (path: string) => {
  const projectName = useProjectStore.getState().currProjectName;
  return await load({ projectName, path, overwrite: false }, `load repository ${path}`);
};

// Loads a tree the user edited in the browser, against the given directory.
export const loadRepositoryJsonIntoProject = async (
  repositoryJSON: object,
  baseDir: string,
  overwrite: boolean,
) => {
  const projectName = useProjectStore.getState().currProjectName;
  return await load({ projectName, repositoryJSON, baseDir, overwrite }, 'add repository');
};

// Loads every registered repository into the current project.
export const loadAllRepositoriesIntoProject = async () => {
  const projectName = useProjectStore.getState().currProjectName;
  return await load({ projectName, overwrite: false }, 'update repositories');
};
