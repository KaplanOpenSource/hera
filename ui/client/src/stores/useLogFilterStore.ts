import { create } from 'zustand';
import { persist } from 'zustand/middleware';
import { LogLineKind } from '../components/workflow/log/classifyLog';

type LogKindVisibility = { [kind in LogLineKind]: boolean };

type LogFilterStore = {
  visible: LogKindVisibility;
  toggle: (kind: LogLineKind) => void;
};

// All kinds start visible except Luigi events, which are hidden by default (they
// are diagnostic noise most of the time; toggle the Events chip to show them).
const defaultVisible = (): LogKindVisibility => {
  return Object.fromEntries(
    Object.values(LogLineKind).map((k) => { return [k, k !== LogLineKind.Event]; })
  ) as LogKindVisibility;
};

// Which workflow-log line kinds are shown. Persisted so the user's choice sticks
// across runs and reloads instead of resetting every time.
export const useLogFilterStore = create<LogFilterStore>()(
  persist(
    (set) => {
      return {
        visible: defaultVisible(),
        toggle: (kind) => {
          return set((state) => { return { visible: { ...state.visible, [kind]: !state.visible[kind] } }; });
        },
      };
    },
    {
      name: 'hera-workflow-log-filters',
      // Merge saved flags over the defaults so a newly added kind stays visible.
      merge: (persisted, current) => {
        const saved = persisted as { visible?: Partial<LogKindVisibility> } | undefined;
        return { ...current, visible: { ...current.visible, ...saved?.visible } };
      },
    }
  )
);
