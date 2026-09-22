import { ProjectEntire } from '@shared/types';
import { create } from 'zustand';
import { ProjectObj } from '../objects/ProjectObj';

export const NO_PROJECT = "* NONE *";
export const DEFAULT_PROJECT = "defaultProject";
export const EMPTY_NAME_PROJECT = "* Empty Name *";

interface ProjectStore {
  currProjectName: string;
  currProject: ProjectEntire | null; // Current project
  selectProject: (newProjectId: string) => void;
  setCurrentProject: (project: ProjectEntire | null) => void; // Sets current project
  getProject: () => ProjectObj | null;
}

export const useProjectStore = create<ProjectStore>((set, get) => ({
  currProjectName: NO_PROJECT,
  currProject: null,
  selectProject: (newProjectName: string) => {
    set({ currProjectName: newProjectName })
  },
  setCurrentProject: (project) => {
    set({ currProject: project })
  },
  getProject: () => {
    const { currProject } = get();
    return currProject ? new ProjectObj(currProject) : null;
  },
}));
