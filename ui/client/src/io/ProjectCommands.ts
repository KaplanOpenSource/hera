import type { PythonCommand } from './fetchPython';
import { DEFAULT_PROJECT } from '../stores/useProjectStore';

export const ProjectCommands = {
  projectNames: (): PythonCommand => ({
    results: ['projects'],
    code: `
from hera.datalayer.project import getProjectList
names = getProjectList()
first = [{"name": n} for n in names if n == "${DEFAULT_PROJECT}"]
rest = [{"name": n} for n in names if n != "${DEFAULT_PROJECT}"]
rest = sorted(rest, key=lambda p: p["name"].lower())
projects = first + rest
`,
  }),
};
