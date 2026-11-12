import { describe, it, expect } from 'vitest';
import { friendlyParamName } from '../src/components/workflow/friendlyParamName';

describe('friendlyParamName', () => {
  it('splits camel case names', () => {
    expect(friendlyParamName('ProjectName')).toBe('Project Name');
    expect(friendlyParamName('projectName')).toBe('Project Name');
  });

  it('splits underscores and dashes', () => {
    expect(friendlyParamName('foo_bar')).toBe('Foo Bar');
    expect(friendlyParamName('foo-bar')).toBe('Foo Bar');
    expect(friendlyParamName('foo__bar')).toBe('Foo Bar');
  });

  it('keeps a run of capitals together', () => {
    expect(friendlyParamName('URL')).toBe('URL');
    expect(friendlyParamName('HTMLPage')).toBe('HTML Page');
  });

  it('capitalizes a single plain word', () => {
    expect(friendlyParamName('name')).toBe('Name');
  });

  it('keeps a name with no words as it is', () => {
    expect(friendlyParamName('')).toBe('');
    expect(friendlyParamName('_')).toBe('_');
  });
});
