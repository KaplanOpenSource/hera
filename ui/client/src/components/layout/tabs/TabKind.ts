// The kinds of tab the dock can hold. The value is also the flexlayout component
// id stored on the tab node, so a node says which class owns it.
export enum TabKind {
  Tree = 'tree',
  Details = 'details',
  Preview = 'preview',
  Canvas = 'canvas',
  Output = 'output',
}
