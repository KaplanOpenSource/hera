# Plan: reference another node's input parameter (#1001, item 1)

Today a node's parameter can only point at another node's **output**. The ticket asks that an input
can be referenced too, so each input row needs a dot on the left (it receives) and one on the right
(it can be read by someone else). Both ways of making the link are in scope: dragging a line, and
the reference completion (right-click menu and inline `{…}` typing).

## What the link must say

An input reference has to be written `{A.Execution.input_parameters.<key>}`. That is the only form
Hermes resolves - a node saves its inputs nested under `Execution`, so `{A.parameters.<key>}` raises
a `KeyError`. No Python change is needed; the dependency scan already treats `{A.…}` as a
requirement on A.

## The work, in five parts

1. **A reference gets a kind.** Output or input. The token writer, the token parser, the edge id and
   the insert/clear helpers all carry it, so an output and a parameter with the same name never mix.

2. **A second dot on each input row.** A source handle at the node's right edge, styled like an
   output's dot, on top-level parameter rows only.

3. **Lines resolve both kinds.** The edge builder matches a referenced name against the other node's
   outputs and its parameter keys, and anchors the line to whichever handle fits. A drag from an
   input's right dot writes an input reference. A node may not reference its own input.

4. **Completion offers inputs too.** The right-click fly-out lists a node's outputs and inputs
   grouped by kind; the inline menu returns options tagged with their kind so a pick writes the
   right section.

5. **Tests and validation.** Cover the new token form, the new edge, the self-reference block and
   the completion options. Then run `ui/client/TEST_UI.md` in order.

## Open question

An input reference reads a value the other node was *given*, not one it produced, yet Hermes still
makes the reader wait for that node to finish. The link is safe but adds ordering the user may not
expect.
