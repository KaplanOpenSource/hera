import { OutputReferenceKind } from './OutputReferenceKind';
import { ReferenceKindRegistry } from './ReferenceKindRegistry';

// The kind every reference is of today.
export const OUTPUT = new OutputReferenceKind();

// The kinds the workflow editor knows how to read and write.
export const knownKinds = new ReferenceKindRegistry([OUTPUT]);
