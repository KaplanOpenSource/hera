import { Sync, SyncDisabled } from '@mui/icons-material';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { DEFAULT_RELOAD_INTERVAL_SECONDS, useViewSettingsStore } from '../../stores/useViewSettingsStore';

// Toolbar toggle for project auto-reload. Off stores null; turning it back on
// restores the default interval.
export const AutoReloadToggle = () => {
  const { viewSettings, setViewSettings } = useViewSettingsStore();
  const autoReloadOff = viewSettings.reloadIntervalSeconds === null;
  return (
    <ButtonTooltip
      title={autoReloadOff
        ? `Turn on auto-reload (${DEFAULT_RELOAD_INTERVAL_SECONDS}s)`
        : 'Turn off auto-reload'}
      onClick={() => setViewSettings({
        reloadIntervalSeconds: autoReloadOff ? DEFAULT_RELOAD_INTERVAL_SECONDS : null,
      })}
      color="inherit"
    >
      {autoReloadOff ? <SyncDisabled /> : <Sync />}
    </ButtonTooltip>
  );
};
