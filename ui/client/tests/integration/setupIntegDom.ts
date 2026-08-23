// jsdom reports every element as 0x0. flexlayout-react lays its tabs out from measured
// rects, so with no size it renders the tab *titles* but never the tab *content* — the
// project tree and the details panel simply do not exist in the DOM. Give jsdom a fixed
// viewport and a ResizeObserver that reports it, so the layout mounts its panels.
const WIDTH = 1280;
const HEIGHT = 900;

Element.prototype.getBoundingClientRect = function (): DOMRect {
  return {
    width: WIDTH, height: HEIGHT,
    top: 0, left: 0, bottom: HEIGHT, right: WIDTH, x: 0, y: 0,
    toJSON() { },
  } as DOMRect;
};

Object.defineProperty(HTMLElement.prototype, 'offsetWidth', { configurable: true, get: () => WIDTH });
Object.defineProperty(HTMLElement.prototype, 'offsetHeight', { configurable: true, get: () => HEIGHT });

class SizedResizeObserver {
  private callback: ResizeObserverCallback;
  constructor(callback: ResizeObserverCallback) { this.callback = callback; }
  observe(target: Element) {
    this.callback(
      [{ target, contentRect: { width: WIDTH, height: HEIGHT } } as ResizeObserverEntry],
      this as unknown as ResizeObserver,
    );
  }
  unobserve() { }
  disconnect() { }
}

global.ResizeObserver = SizedResizeObserver as unknown as typeof ResizeObserver;
