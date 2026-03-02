import { describe, it, expect } from 'vitest';
import { resolveProjectFromUrl } from '../src/io/FetchProjects';

const projects = [
  { name: 'Alpha' },
  { name: 'Beta' },
  { name: 'My Project' },
];

describe('resolveProjectFromUrl', () => {
  it('returns matching project name for a valid encoded URL param', () => {
    expect(resolveProjectFromUrl(encodeURIComponent('My Project'), projects)).toBe('My Project');
  });

  it('returns undefined when URL param does not match any project', () => {
    expect(resolveProjectFromUrl('NonExistent', projects)).toBeUndefined();
  });
});
