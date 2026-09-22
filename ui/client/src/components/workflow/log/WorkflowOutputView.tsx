import { Box, CircularProgress, Typography } from '@mui/material';
import { useCallback, useEffect, useRef, useState } from 'react';
import { chunksToText, WorkflowChunk } from '../../../io/runWorkflow';
import { useLogFilterStore } from '../../../stores/useLogFilterStore';
import { useWorkflowFocusStore } from '../../../stores/useWorkflowFocusStore';
import { chunkedMetrics, flatMetrics, LogMetrics } from './logMetrics';
import { LogToolbar } from './LogToolbar';
import { isTaskChunk } from './WorkflowChunkLog';
import { WorkflowChunkedLog } from './WorkflowChunkedLog';
import { WorkflowLogView } from './WorkflowLogView';
import { WorkflowNodeSelect } from './WorkflowNodeSelect';
import { ChunkEdge, currentChunkIndex } from './visibleChunk';

// Shows a workflow run's output as per-task cards, growing live as the run streams
// (with a small "running" hint) and staying in the same shape once it finishes, so
// the view does not reflow at the end. On failure it also shows the error. The
// log-level filter + copy-all toolbar sits in a bottom bar, outside the scrolling
// log, so it stays put while the log scrolls.
export const WorkflowOutputView = ({
  running,
  chunks,
  error,
  workflowName,
}: {
  running: boolean,
  chunks?: WorkflowChunk[] | null,
  error: string | null,
  workflowName?: string,
}) => {
  const visible = useLogFilterStore((state) => { return state.visible; });
  const toggle = useLogFilterStore((state) => { return state.toggle; });
  const focusNode = useWorkflowFocusStore((state) => { return state.focusNode; });

  const chunkList = chunks ?? [];
  // Cards throughout: the chunks carry their task name while running too, so the
  // view grows live and does not switch shape when the run finishes.
  const showChunked = chunkList.length > 0;
  // Show the log while running and whenever there is output, even on failure: the
  // log that led to the error stays visible (below the error message). Only when a
  // failure produced no output at all is there nothing to show but the error.
  const showLog = running || chunkList.length > 0;

  const scrollRef = useRef<HTMLDivElement>(null);
  const [currentIndex, setCurrentIndex] = useState<number | undefined>(undefined);

  // The rendered cards, in document order, marked by WorkflowChunkLog.
  const cards = (): HTMLElement[] => {
    return Array.from(scrollRef.current?.querySelectorAll<HTMLElement>('[data-chunk-index]') ?? []);
  };

  const syncCurrent = useCallback(() => {
    const view = scrollRef.current;
    if (!view) {
      return;
    }
    const viewTop = view.getBoundingClientRect().top;
    // Only node cards; the select cannot show the setup / between / final filler.
    const edges: ChunkEdge[] = cards()
      .map((card) => {
        const rect = card.getBoundingClientRect();
        return { index: Number(card.dataset.chunkIndex), top: rect.top, bottom: rect.bottom };
      })
      .filter((edge) => { return isTaskChunk(chunkList[edge.index]?.name ?? ''); });
    setCurrentIndex(currentChunkIndex(viewTop, edges));
  }, [chunkList]);

  // A live run adds cards, so the shown node is recomputed as the log grows too.
  useEffect(syncCurrent, [syncCurrent, chunkList.length]);

  const scrollToChunk = (index: number) => {
    const view = scrollRef.current;
    const card = cards().find((c) => { return Number(c.dataset.chunkIndex) === index; });
    if (!view || !card) {
      return;
    }
    view.scrollTop += card.getBoundingClientRect().top - view.getBoundingClientRect().top;
    setCurrentIndex(index);
    // The select only offers task chunks, so the name is a node on the canvas.
    if (workflowName) {
      focusNode(workflowName, chunkList[index].name);
    }
  };

  let metrics: LogMetrics;
  if (showChunked) {
    metrics = chunkedMetrics(chunkList);
  } else {
    metrics = flatMetrics(chunksToText(chunkList));
  }

  return (
    <Box sx={{ position: 'relative', height: '100%', display: 'flex', flexDirection: 'column', minHeight: 0 }}>
      {/* Floats over the log rather than taking a bar of its own. */}
      {showChunked && (
        <Box sx={{ position: 'absolute', top: 8, right: 16, zIndex: 2 }}>
          <WorkflowNodeSelect chunks={chunkList} currentIndex={currentIndex} onPick={scrollToChunk} />
        </Box>
      )}
      <Box ref={scrollRef} onScroll={syncCurrent} sx={{ flexGrow: 1, overflow: 'auto', p: 1, minHeight: 0 }}>
        {!running && error && (
          <Typography color="error" sx={{ whiteSpace: 'pre-wrap', wordBreak: 'break-word' }}>
            {error}
          </Typography>
        )}
        {showLog && (
          <>
            {running && (
              <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, pb: 1 }}>
                <CircularProgress size={16} />
                <Typography variant="body2" color="text.secondary">Running…</Typography>
              </Box>
            )}
            {showChunked
              ? <WorkflowChunkedLog chunks={chunkList} />
              : <WorkflowLogView output={metrics.fullText} />}
          </>
        )}
      </Box>
      {showLog && (
        <Box sx={{ display: 'flex', borderTop: 1, borderColor: 'divider', px: 1, py: 0.5 }}>
          <LogToolbar counts={metrics.counts} visible={visible} onToggle={toggle} fullText={metrics.fullText} />
        </Box>
      )}
    </Box>
  );
};
