import { describe, it, expect } from 'vitest';
import { CanvasResize, canvasResizeAction } from '../src/components/workflow/canvasResize';

describe('canvasResizeAction', () => {
  it('does nothing on the first observation', () => {
    expect(canvasResizeAction(null, { width: 800, height: 600 })).toBe(CanvasResize.None);
  });

  it('does nothing when a dimension is zero', () => {
    expect(canvasResizeAction({ width: 800, height: 600 }, { width: 0, height: 600 })).toBe(CanvasResize.None);
    expect(canvasResizeAction({ width: 0, height: 600 }, { width: 800, height: 600 })).toBe(CanvasResize.None);
  });

  it('does nothing when the size did not change', () => {
    expect(canvasResizeAction({ width: 800, height: 600 }, { width: 800, height: 600 })).toBe(CanvasResize.None);
  });

  it('fits when a tab docks beside the canvas and narrows it', () => {
    expect(canvasResizeAction({ width: 800, height: 600 }, { width: 500, height: 600 })).toBe(CanvasResize.Fit);
  });

  it('fits when the canvas widens again', () => {
    expect(canvasResizeAction({ width: 500, height: 600 }, { width: 800, height: 600 })).toBe(CanvasResize.Fit);
  });

  it('fits on maximize, where both sides change', () => {
    expect(canvasResizeAction({ width: 500, height: 300 }, { width: 1200, height: 900 })).toBe(CanvasResize.Fit);
  });

  it('scales the zoom when only the height changes', () => {
    expect(canvasResizeAction({ width: 800, height: 600 }, { width: 800, height: 300 })).toBe(CanvasResize.ScaleZoom);
  });
});
