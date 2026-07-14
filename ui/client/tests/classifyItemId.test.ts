import { describe, it, expect } from 'vitest';
import {
  CENTRAL_REPO_FOLDER_ID,
  classifyItemId,
  idDocId,
  idRepoId,
  ItemKind,
} from '../src/shared/idDocId';

describe('classifyItemId', () => {
  it('classifies the central repo folder id', () => {
    expect(classifyItemId(CENTRAL_REPO_FOLDER_ID)).toBe(ItemKind.CentralRepo);
  });

  it('classifies split ids', () => {
    expect(classifyItemId('split_/toolkit=MyToolkit')).toBe(ItemKind.Split);
    expect(classifyItemId('split_field=value')).toBe(ItemKind.Split);
  });

  it('classifies document ids', () => {
    expect(classifyItemId(idDocId('abc123'))).toBe(ItemKind.Document);
  });

  it('classifies repo ids', () => {
    expect(classifyItemId(idRepoId('my_repo'))).toBe(ItemKind.Repo);
  });

  it('falls back to config for unrecognized ids', () => {
    expect(classifyItemId('config')).toBe(ItemKind.Config);
    expect(classifyItemId('anything-else')).toBe(ItemKind.Config);
  });

  it('checks split before document and repo (prefixes are distinct)', () => {
    // A split id never looks like a document or repo id, and vice versa.
    expect(classifyItemId('split_/toolkit=X')).toBe(ItemKind.Split);
    expect(classifyItemId(idDocId('split_looking'))).toBe(ItemKind.Document);
  });
});
