import { Switch, Tooltip } from '@mui/material';
import { DEFAULT_RELOAD_INTERVAL_SECONDS, useViewSettingsStore } from '../../stores/useViewSettingsStore';

// Toolbar toggle switch for project auto-reload. Off stores null; turning it back
// on restores the default interval. The "on" color is the cyan accent.
export const AutoReloadToggle = () => {
  const { viewSettings, setViewSettings } = useViewSettingsStore();
  const autoReloadOn = viewSettings.reloadIntervalSeconds !== null;
  return (
    <Tooltip
      title={autoReloadOn
        ? 'Auto-reload is on, click to turn it off'
        : `Auto-reload is off, click to set it on to ${DEFAULT_RELOAD_INTERVAL_SECONDS}s`}
    >
      <Switch
        checked={autoReloadOn}
        onChange={() => setViewSettings({
          reloadIntervalSeconds: autoReloadOn ? null : DEFAULT_RELOAD_INTERVAL_SECONDS,
        })}
        sx={{
          '& .MuiSwitch-switchBase.Mui-checked': { color: '#22d3ee' },
          '& .MuiSwitch-switchBase.Mui-checked + .MuiSwitch-track': { backgroundColor: '#22d3ee' },
        }}
      />
    </Tooltip>
  );
};
