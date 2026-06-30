import { ProjectEntire, ProjectName, Toolkit } from '@shared/types';
import { create } from 'zustand';
import { ProjectObj } from '../objects/ProjectObj';
import { ToolkitObj } from '../objects/ToolkitObj';

export const NO_PROJECT = "* NONE *";
export const DEFAULT_PROJECT = "defaultProject";
export const EMPTY_NAME_PROJECT = "* Empty Name *";

interface ProjectStore {
  projectNames: ProjectName[]; // List of project names
  currProjectName: string;
  currProject: ProjectEntire | null; // Current project
  toolkits: Toolkit[];
  setProjectNames: (names: ProjectName[]) => void; // Sets project names
  selectProject: (newProjectId: string) => void;
  setCurrentProject: (project: ProjectEntire | null) => void; // Sets current project
  setToolkits: (val: Toolkit[]) => void;
  getProject: () => ProjectObj | null;
  getProjectToolkitKeys: () => string[];
}

export const useProjectStore = create<ProjectStore>((set, get) => ({
  projectNames: [],
  currProjectName: NO_PROJECT,
  currProject: null,
  toolkits: [],
  setProjectNames: (names) => {
    set({ projectNames: names })
  },
  selectProject: (newProjectName: string) => {
    set({ currProjectName: newProjectName })
  },
  setCurrentProject: (project) => {
    set({ currProject: project })
  },
  setToolkits: (val) => {
    set({ toolkits: val })
  },
  getProject: () => {
    const { currProject } = get();
    return currProject ? new ProjectObj(currProject) : null;
  },
  getProjectToolkitKeys: () => {
    const { toolkits } = get();
    const documents = get().getProject()?.documents ?? [];
    const docToolkitNames = [...new Set(
      documents.map(d => d.toolkit).filter(Boolean) as string[]
    )];
    return toolkits
      .filter(t => docToolkitNames.some(dt =>
        new ToolkitObj(t).matches(dt)
      ))
      .map(t => t.toolkit);
  },
}));
