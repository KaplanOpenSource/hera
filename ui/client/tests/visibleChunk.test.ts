import { describe, it, expect } from 'vitest';
import { currentChunkIndex } from '../src/components/workflow/log/visibleChunk';

// Only node cards reach here; the filler between them leaves the gaps.
// Node 0 spans 0-80, node 1 spans 200-300, node 2 spans 420-500.
const edges = [
  { index: 0, top: 0, bottom: 80 },
  { index: 1, top: 200, bottom: 300 },
  { index: 2, top: 420, bottom: 500 },
];

describe('currentChunkIndex', () => {
  it('has nothing to show with no node cards', () => {
    expect(currentChunkIndex(100, [])).toBeUndefined();
  });

  it('picks the node the view sits on', () => {
    expect(currentChunkIndex(0, edges)).toBe(0);
    expect(currentChunkIndex(250, edges)).toBe(1);
  });

  it('picks the upper node when two span the top edge', () => {
    const overlapping = [
      { index: 0, top: 0, bottom: 300 },
      { index: 1, top: 100, bottom: 400 },
    ];
    expect(currentChunkIndex(150, overlapping)).toBe(0);
  });

  it('picks the nearest node when the view sits on filler above it', () => {
    // 180 is in the gap after node 0: 100 past its end, 20 before node 1.
    expect(currentChunkIndex(180, edges)).toBe(1);
  });

  it('picks the nearest node when the view sits on filler below it', () => {
    // 110 is in the gap after node 0: 30 past its end, 90 before node 1.
    expect(currentChunkIndex(110, edges)).toBe(0);
  });

  it('keeps the last node once everything is scrolled past', () => {
    expect(currentChunkIndex(900, edges)).toBe(2);
  });

  it('treats a node ending exactly at the top edge as scrolled past', () => {
    // 80 is node 0's bottom, so it is behind; node 1 starts 120 below.
    expect(currentChunkIndex(80, edges)).toBe(0);
  });
});
