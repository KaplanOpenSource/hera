import { LogLineKind } from './classifyLog';

// Shared per-kind label and text color, used by both the log lines and the
// filter buttons so they always match. Colors are theme palette keys.
export interface LogKindStyle {
  label: string,
  color: string,
  bold?: boolean,
}

export const KIND_STYLE: { [kind in LogLineKind]: LogKindStyle } = {
  [LogLineKind.Debug]: { label: 'Debug', color: 'text.disabled' },
  [LogLineKind.Info]: { label: 'Info', color: 'primary.main' },
  [LogLineKind.Warning]: { label: 'Warning', color: 'warning.main' },
  [LogLineKind.Error]: { label: 'Error', color: 'error.main', bold: true },
  [LogLineKind.Summary]: { label: 'Summary', color: 'secondary.main', bold: true },
  [LogLineKind.Output]: { label: 'Output', color: 'success.main', bold: true },
};

// Display order for the filter buttons.
export const KIND_ORDER: LogLineKind[] = [
  LogLineKind.Output,
  LogLineKind.Summary,
  LogLineKind.Info,
  LogLineKind.Debug,
  LogLineKind.Warning,
  LogLineKind.Error,
];
