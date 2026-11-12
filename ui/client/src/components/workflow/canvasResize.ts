// What the canvas should do to its viewport when its size changes.
export enum CanvasResize {
  // Size is not usable yet (first observation, or a zero dimension).
  None = 'none',
  // Fit the whole graph: the width changed, so narrowing would cut nodes off.
  Fit = 'fit',
  // Height only: scale the zoom by the same ratio to keep the same slice framed.
  ScaleZoom = 'scaleZoom',
}

export type CanvasSize = { width: number, height: number };

export const canvasResizeAction = (prev: CanvasSize | null, next: CanvasSize): CanvasResize => {
  if (!prev || !next.width || !next.height || !prev.width || !prev.height) {
    return CanvasResize.None;
  }
  if (prev.width !== next.width) {
    return CanvasResize.Fit;
  }
  if (prev.height !== next.height) {
    return CanvasResize.ScaleZoom;
  }
  return CanvasResize.None;
};
