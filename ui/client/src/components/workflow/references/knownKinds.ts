import { InputReferenceKind } from './InputReferenceKind';
import { OutputReferenceKind } from './OutputReferenceKind';
import { ReferenceKindRegistry } from './ReferenceKindRegistry';

// A reference to a node's output.
export const OUTPUT = new OutputReferenceKind();

// A reference to a node's own input parameter.
export const INPUT = new InputReferenceKind();

// The kinds the workflow editor knows how to read and write.
export const knownKinds = new ReferenceKindRegistry([OUTPUT, INPUT]);
