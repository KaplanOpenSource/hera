import { useEffect } from 'react';
import { NO_PROJECT, useProjectStore } from '../stores/useProjectStore';
import { useViewSettingsStore } from '../stores/useViewSettingsStore';
import { fetchProjectDetails } from './FetchProjects';

// Periodically reloads the open project's documents into the store (silently),
// so the tree and open documents — which read from the store — stay up to date.
// The interval is configurable in settings; 0 (or null) turns auto-reload off.
export const ProjectAutoReload = () => {
  const currProjectName = useProjectStore(s => s.currProjectName);
  const reloadIntervalSeconds = useViewSettingsStore(s => s.viewSettings.reloadIntervalSeconds);

  useEffect(() => {
    if (currProjectName === NO_PROJECT || !reloadIntervalSeconds) {
      return;
    }
    const interval = setInterval(() => {
      fetchProjectDetails(currProjectName, true);
    }, reloadIntervalSeconds * 1000);
    return () => {
      clearInterval(interval);
    };
  }, [currProjectName, reloadIntervalSeconds]);

  return null;
};
