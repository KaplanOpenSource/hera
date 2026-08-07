// Classifies raw workflow console output into indexed, typed lines so the UI can
// color each one. Kept entirely on the client — the server returns the raw log.

export enum LogLineKind {
  Debug = 'DEBUG',
  Info = 'INFO',
  Warning = 'WARNING',
  Error = 'ERROR',
  Summary = 'SUMMARY',
  Output = 'OUTPUT',
}

export interface ClassifiedLine {
  index: number,
  kind: LogLineKind,
  text: string,
}

// Leading "LEVEL:" prefixes emitted by Luigi's logger.
const PREFIX_TO_KIND: { [prefix: string]: LogLineKind } = {
  'DEBUG:': LogLineKind.Debug,
  'INFO:': LogLineKind.Info,
  'WARNING:': LogLineKind.Warning,
  'ERROR:': LogLineKind.Error,
  'CRITICAL:': LogLineKind.Error,
};

const SUMMARY_MARKER = 'Luigi Execution Summary';

const kindForLine = (
  text: string,
): LogLineKind | null => {
  const prefix = Object.keys(PREFIX_TO_KIND).find((p) => { return text.startsWith(p); });
  return prefix ? PREFIX_TO_KIND[prefix] : null;
};

// Lines between the two "===== Luigi Execution Summary =====" markers (and the
// markers themselves) are the summary block. Anything without a level prefix and
// outside that block is the task's own stdout (Output).
export const classifyLog = (
  raw: string,
): ClassifiedLine[] => {
  let insideSummary = false;
  return raw.split('\n').map((text, index) => {
    const isMarker = text.includes(SUMMARY_MARKER);
    if (isMarker) {
      insideSummary = !insideSummary;
      return { index, kind: LogLineKind.Summary, text };
    }
    if (insideSummary) {
      return { index, kind: LogLineKind.Summary, text };
    }
    const prefixKind = kindForLine(text);
    return { index, kind: prefixKind ?? LogLineKind.Output, text };
  });
};
