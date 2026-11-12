# Plan: reference another node's input parameter (#1001, item 1)

A node's parameter can only point at another node's output today. It should be able to point at
another node's input too. So an input row gets a dot on both sides: the left one takes a value in,
the right one lets another node read it. Dragging a line and the completion menus both have to work.

An input reference is written `{A.Execution.input_parameters.<key>}`. That is the only form Hermes
resolves. No Python change.

## Already there

`src/components/workflow/references/` holds one class per reference kind, and the registry reads a
written token and says which kind it is. Input rows have their left dot; output rows have their
right dot.

## The work

1. Add the input kind. It cannot use the `in` handle mark, which the left dot already owns.
2. Put a right dot on each top-level input row. This needs a new render slot after the value editor.
3. Stop a node pointing at its own input.
4. Teach the inline `{…}` typing which section is being typed. It only tracks the node name today.
5. List both kinds in the right-click menu.
6. Tests, then `ui/client/TEST_UI.md` in order.

Step 2 is the only visible change. Step 4 is the real work.

## Note

An input reference makes the reader wait for the other node to finish, even though the value was
known before the run. That is Hermes' run order, not the client's or the server's. Nothing to do
here.
