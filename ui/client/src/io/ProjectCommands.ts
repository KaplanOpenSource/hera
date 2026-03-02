import type { PythonCommand } from './fetchPython';

export const ProjectCommands = {
  projectNames: (): PythonCommand => ({
    results: ['projects'],
    code: `
from hera.datalayer.project import getProjectList
names = getProjectList()
projects = [{"name": proj} for proj in names]
`,
  }),
};
