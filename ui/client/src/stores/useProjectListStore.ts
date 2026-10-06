import { ProjectName } from '@shared/types';
import { create } from 'zustand';

interface ProjectListStore {
  projectNames: ProjectName[];
  setProjectNames: (names: ProjectName[]) => void;
}

export const useProjectListStore = create<ProjectListStore>((set) => ({
  projectNames: [],
  setProjectNames: (names) => {
    set({ projectNames: names })
  },
}));
