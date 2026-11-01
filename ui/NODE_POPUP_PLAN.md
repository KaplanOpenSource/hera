# Node collapse rework (based on #1134)

## Goal

- Nodes always stay small. Hovering does not open them.
- The pin button opens a node in place, as today. It pushes the other nodes aside.
- A new button opens the node in a pop-up, where it can be edited.

## Steps

1. Stop opening nodes on hover.
2. Show the node's buttons (pin, pop-up, delete) on the small card too.
3. Make the node editor usable outside the canvas.
4. Add the pop-up button and the pop-up window.
5. Add tests and run the UI checklist.

## Limit

- Lines between nodes can't be dragged from the pop-up. Use a pinned node for that.
