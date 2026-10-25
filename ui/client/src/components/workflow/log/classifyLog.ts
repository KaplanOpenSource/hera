// Classifies raw workflow console output into indexed, typed lines so the UI can
// color each one. Kept entirely on the client — the server returns the raw log.

export enum LogLineKind {
  Debug = 'DEBUG',
  Info = 'INFO',
  Warning = 'WARNING',
  Error = 'ERROR',
  Summary = 'SUMMARY',
  Output = 'OUTPUT',
  Technical = 'TECHNICAL',
  Event = 'EVENT',
}

// Prefix every Luigi task/dependency event line carries. Emitted by the server's
// in-process runner (execute_workflow_inprocess.py, EVENT_PREFIX). These are hidden
// by default so they don't clutter the task's own output.
export const EVENT_PREFIX = '[luigi-event]';

export interface ClassifiedLine {
  index: number,
  kind: LogLineKind,
  text: string,
}

// Leading "LEVEL:" prefixes emitted by Python/Luigi's logger. CRITICAL is treated
// as Error and WARN as Warning (they're synonyms).
const PREFIX_TO_KIND: { [prefix: string]: LogLineKind } = {
  'DEBUG:': LogLineKind.Debug,
  'INFO:': LogLineKind.Info,
  'WARNING:': LogLineKind.Warning,
  'WARN:': LogLineKind.Warning,
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
// markers themselves) are the summary block. Once that block closes, everything
// after it is Technical (Hermes loading node outputs, load problems, etc.).
// Before the block, a line without a level prefix is the task's own stdout (Output).
export const classifyLog = (
  raw: string,
): ClassifiedLine[] => {
  let insideSummary = false;
  let afterSummary = false;
  return raw.split('\n').map((text, index) => {
    // Luigi event lines are self-identifying by prefix; classify them first so the
    // summary/output logic below never reclassifies them.
    if (text.startsWith(EVENT_PREFIX)) {
      return { index, kind: LogLineKind.Event, text };
    }
    const isMarker = text.includes(SUMMARY_MARKER);
    if (isMarker) {
      if (insideSummary) {
        afterSummary = true;
      }
      insideSummary = !insideSummary;
      return { index, kind: LogLineKind.Summary, text };
    }
    if (insideSummary) {
      return { index, kind: LogLineKind.Summary, text };
    }
    if (afterSummary) {
      return { index, kind: LogLineKind.Technical, text };
    }
    const prefixKind = kindForLine(text);
    return { index, kind: prefixKind ?? LogLineKind.Output, text };
  });
};
