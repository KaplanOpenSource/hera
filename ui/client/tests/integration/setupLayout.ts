// flexlayout-react renders its tab contents only after it measures a non-zero
// container. jsdom reports every element as 0x0 and never fires ResizeObserver,
// so the project tree never mounts. Give the layout a real viewport.
const WIDTH = 1280;
const HEIGHT = 800;

const rect: DOMRect = {
  width: WIDTH, height: HEIGHT,
  top: 0, left: 0, bottom: HEIGHT, right: WIDTH, x: 0, y: 0,
  toJSON: () => ({}),
};

Element.prototype.getBoundingClientRect = () => rect;

Object.defineProperty(HTMLElement.prototype, 'offsetWidth', { configurable: true, value: WIDTH });
Object.defineProperty(HTMLElement.prototype, 'offsetHeight', { configurable: true, value: HEIGHT });
Object.defineProperty(HTMLElement.prototype, 'clientWidth', { configurable: true, value: WIDTH });
Object.defineProperty(HTMLElement.prototype, 'clientHeight', { configurable: true, value: HEIGHT });

class SizedResizeObserver {
  private callback: ResizeObserverCallback;

  constructor(callback: ResizeObserverCallback) {
    this.callback = callback;
  }

  observe(target: Element) {
    this.callback([{ target, contentRect: rect } as ResizeObserverEntry], this as unknown as ResizeObserver);
  }

  unobserve() {}
  disconnect() {}
}

global.ResizeObserver = SizedResizeObserver as unknown as typeof ResizeObserver;
