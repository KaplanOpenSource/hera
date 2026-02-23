import { create } from 'zustand';

export type ViewSettingsType = {
  minGroupSize: number;
  maxDepth: number;
  maxBranches: number;
  showEmptyToolkits: boolean;
  showDocumentPreview: boolean;
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
  showEmptyToolkits: false,
  showDocumentPreview: true,
};

export const useViewSettingsStore = create<ViewSettingsStore>((set) => ({
  viewSettings: defaultSettings,
  setViewSettings: (settings) =>
    set((state) => ({ viewSettings: { ...state.viewSettings, ...settings } })),
  resetViewSettings: () => set({ viewSettings: defaultSettings }),
}));
