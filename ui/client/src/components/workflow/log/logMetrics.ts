import { WorkflowChunk } from '../../../io/runWorkflow';
import { classifyLog, EVENT_PREFIX, LogLineKind } from './classifyLog';

// Drops the [luigi-event] marker lines: the task name is shown as a header/label,
// so the raw event markers are noise in the displayed log.
export const withoutEventLines = (
  text: string,
): string => {
  return text
    .split('\n')
    .filter((line) => { return !line.startsWith(EVENT_PREFIX); })
    .join('\n');
};

const emptyCounts = (): { [kind in LogLineKind]: number } => {
  return Object.fromEntries(Object.values(LogLineKind).map((k) => { return [k, 0]; })) as {
    [kind in LogLineKind]: number
  };
};

// Line counts per kind and the plain text for copy-all, feeding the shared toolbar.
export type LogMetrics = {
  counts: { [kind in LogLineKind]: number },
  fullText: string,
};

// Metrics for the grouped per-task view: event markers stripped, and only chunks
// with visible content contribute (empty ones render no card).
export const chunkedMetrics = (
  chunks: WorkflowChunk[],
): LogMetrics => {
  const counts = emptyCounts();
  chunks.forEach((chunk) => {
    const lines = classifyLog(withoutEventLines(chunk.text));
    if (lines.some((line) => { return line.text.trim() !== ''; })) {
      lines.forEach((line) => { counts[line.kind] += 1; });
    }
  });
  const fullText = chunks.map((chunk) => { return withoutEventLines(chunk.text); }).join('');
  return { counts, fullText };
};

// Metrics for the flat running view: every line of the raw text counts.
export const flatMetrics = (
  text: string,
): LogMetrics => {
  const counts = emptyCounts();
  classifyLog(text).forEach((line) => { counts[line.kind] += 1; });
  return { counts, fullText: text };
};
