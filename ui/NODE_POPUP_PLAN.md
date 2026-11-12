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

## Done

- Steps 3-5 landed in b229eb1c9, merged to master in #1138.
- Steps 1-2 did not land. Hover still opens the editor. We keep it that way.

## Next (new info on #1134, Nov 4)

- Lior is fine with a small node about 50% taller. Try it.

1. Show the outputs on the small card. Done.
2. Put a thin divider between the inputs and the outputs. Done.
3. Give the small card the same line dots as the editor, so lines reach the rows. Done.
   A line into a nested value lands on its top-level row. The dots hide while hover opens the editor.
4. Add tests, run the UI checklist, post a screenshot on #1134.
