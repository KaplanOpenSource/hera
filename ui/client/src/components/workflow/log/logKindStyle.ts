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
  [LogLineKind.Technical]: { label: 'Technical', color: 'info.main' },
};

// Display order for the filter buttons: task output first, then the log levels by
// severity (high to low, so Debug is the rightmost level), then Summary and
// Technical on the right.
export const KIND_ORDER: LogLineKind[] = [
  LogLineKind.Output,
  LogLineKind.Error,
  LogLineKind.Warning,
  LogLineKind.Info,
  LogLineKind.Debug,
  LogLineKind.Summary,
  LogLineKind.Technical,
];
