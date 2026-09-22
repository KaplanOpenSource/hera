// A rendered card: its chunk index and where its edges sit on screen.
export type ChunkEdge = { index: number, top: number, bottom: number };

// The node card to mark as current: the one at the top of the scroll view, or,
// when the view sits on filler between nodes, the nearest node card to it. With
// two node cards spanning the top edge this picks the upper one.
export const currentChunkIndex = (viewTop: number, edges: ChunkEdge[]): number | undefined => {
  let nearest: ChunkEdge | undefined;
  let nearestGap = Infinity;
  for (const edge of edges) {
    if (edge.top <= viewTop && viewTop < edge.bottom) {
      return edge.index;
    }
    const gap = edge.top > viewTop ? edge.top - viewTop : viewTop - edge.bottom;
    if (gap < nearestGap) {
      nearest = edge;
      nearestGap = gap;
    }
  }
  return nearest?.index;
};
