import { ProjectEntire, ProjectName } from '@shared/types';
import { create } from 'zustand';

export const NO_PROJECT = "* NONE *";

interface ProjectStore {
  projectNames: ProjectName[]; // List of project names
  currProjectName: string;
  currProject: ProjectEntire | null; // Current project
  setProjectNames: (names: ProjectName[]) => void; // Sets project names
  selectProject: (newProjectId: string) => void;
  setCurrentProject: (project: ProjectEntire | null) => void; // Sets current project
}

export const useProjectStore = create<ProjectStore>((set) => ({
  projectNames: [],
  currProjectName: NO_PROJECT,
  currProject: null,
  setProjectNames: (names) => set({ projectNames: names }),
  selectProject: (newProjectName: string) => set({ currProjectName: newProjectName }),
  setCurrentProject: (project) => set({ currProject: project }),
}));
