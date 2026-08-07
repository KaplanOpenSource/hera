import { create } from 'zustand';
import { persist } from 'zustand/middleware';

export enum ThemeMode {
  Light = 'light',
  Dark = 'dark',
}

export type ViewSettingsType = {
  minGroupSize: number;
  maxDepth: number;
  maxBranches: number;
  firstBranchHeadFields: boolean;
  showDocumentPreview: boolean;
  // How often (seconds) to auto-reload the open project. null turns auto-reload off.
  reloadIntervalSeconds: number | null;
  themeMode: ThemeMode;
};

type ViewSettingsStore = {
  viewSettings: ViewSettingsType;
  setViewSettings: (settings: Partial<ViewSettingsType>) => void;
  resetViewSettings: () => void;
};

// Default auto-reload interval (seconds) used when auto-reload is turned back on.
export const DEFAULT_RELOAD_INTERVAL_SECONDS = 5;

const defaultSettings: ViewSettingsType = {
  minGroupSize: 2,
  maxDepth: 5,
  maxBranches: 50,
  firstBranchHeadFields: true,
  showDocumentPreview: true,
  reloadIntervalSeconds: DEFAULT_RELOAD_INTERVAL_SECONDS,
  themeMode: ThemeMode.Light,
};

export const useViewSettingsStore = create<ViewSettingsStore>()(
  persist(
    (set) => ({
      viewSettings: defaultSettings,
      setViewSettings: (settings) =>
        set((state) => ({ viewSettings: { ...state.viewSettings, ...settings } })),
      resetViewSettings: () =>
        set((state) => ({
          // Keep the chosen theme; reset only the tree/view knobs.
          viewSettings: { ...defaultSettings, themeMode: state.viewSettings.themeMode },
        })),
    }),
    {
      name: 'hera-view-settings',
      // Persist all view settings (everything shown in the Settings dialog) so
      // the user's choices survive reloads. Only the store functions are dropped.
      partialize: (state) => ({ viewSettings: state.viewSettings }),
      merge: (persisted, current) => {
        const saved = persisted as { viewSettings?: Partial<ViewSettingsType> } | undefined;
        return {
          ...current,
          viewSettings: { ...current.viewSettings, ...saved?.viewSettings },
        };
      },
    }
  )
);
