import { create } from 'zustand';

export type ViewSettingsType = {
  minGroupSize: number;
  maxDepth: number;
  maxBranches: number;
  firstBranchHeadFields: boolean;
  showDocumentPreview: boolean;
  // How often (seconds) to auto-reload the open project. null turns auto-reload off.
  reloadIntervalSeconds: number | null;
};

type ViewSettingsStore = {
  viewSettings: ViewSettingsType;
  setViewSettings: (settings: Partial<ViewSettingsType>) => void;
  resetViewSettings: () => void;
};

const defaultSettings: ViewSettingsType = {
  minGroupSize: 2,
  maxDepth: 5,
  maxBranches: 50,
  firstBranchHeadFields: true,
  showDocumentPreview: true,
  reloadIntervalSeconds: 5,
};

export const useViewSettingsStore = create<ViewSettingsStore>((set) => ({
  viewSettings: defaultSettings,
  setViewSettings: (settings) =>
    set((state) => ({ viewSettings: { ...state.viewSettings, ...settings } })),
  resetViewSettings: () => set({ viewSettings: defaultSettings }),
}));
