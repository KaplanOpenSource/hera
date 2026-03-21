import { describe, it, expect } from 'vitest';
import { ProjectObj, DocumentObj } from '../src/objects/ProjectObj';
import { ProjectEntire } from '../src/shared/types';

const makeProject = (name: string, docs: any[] = []): ProjectEntire => ({
  name,
  documents: docs,
});

const makeDoc = (overrides: Partial<any> = {}) => ({
  _cls: 'Metadata.Cache',
  _id: { $oid: 'abc123' },
  projectName: 'TestProject',
  desc: { datasourceName: 'myDoc' },
  type: 'SomeType',
  resource: '',
  dataFormat: 'string',
  ...overrides,
});

describe('ProjectObj', () => {
  it('returns the project name', () => {
    const p = new ProjectObj(makeProject('Alpha'));
    expect(p.name).toBe('Alpha');
  });

  it('allDocuments includes config documents', () => {
    const docs = [
      makeDoc({ type: 'TestProject__config__' }),
      makeDoc({ _id: { $oid: 'def456' }, type: 'RegularDoc' }),
    ];
    const p = new ProjectObj(makeProject('TestProject', docs));
    expect(p.allDocuments).toHaveLength(2);
  });

  it('documents filters out config documents', () => {
    const docs = [
      makeDoc({ type: 'TestProject__config__' }),
      makeDoc({ _id: { $oid: 'def456' }, type: 'RegularDoc' }),
    ];
    const p = new ProjectObj(makeProject('TestProject', docs));
    expect(p.documents).toHaveLength(1);
    expect(p.documents[0].data.type).toBe('RegularDoc');
  });

  it('configDocument returns the config doc', () => {
    const docs = [
      makeDoc({ type: 'TestProject__config__' }),
      makeDoc({ _id: { $oid: 'def456' }, type: 'RegularDoc' }),
    ];
    const p = new ProjectObj(makeProject('TestProject', docs));
    expect(p.configDocument).toBeDefined();
    expect(p.configDocument!.isConfig).toBe(true);
  });

  it('configDocument is undefined when no config doc exists', () => {
    const docs = [makeDoc({ type: 'RegularDoc' })];
    const p = new ProjectObj(makeProject('TestProject', docs));
    expect(p.configDocument).toBeUndefined();
  });
});

describe('DocumentObj', () => {
  const project = new ProjectObj(makeProject('TestProject', []));

  it('docid returns the OID', () => {
    const d = new DocumentObj(makeDoc({ _id: { $oid: 'xyz789' } }), project);
    expect(d.docid).toBe('xyz789');
  });

  it('name prefers datasourceName', () => {
    const d = new DocumentObj(makeDoc({ desc: { datasourceName: 'MyDS' }, type: 'T' }), project);
    expect(d.name).toBe('MyDS');
  });

  it('name falls back to type when no datasourceName', () => {
    const d = new DocumentObj(makeDoc({ desc: {}, type: 'FallbackType' }), project);
    expect(d.name).toBe('FallbackType');
  });

  it('name falls back to _cls when no datasourceName or type', () => {
    const d = new DocumentObj(makeDoc({ desc: {}, type: '', _cls: 'Meta.X' }), project);
    expect(d.name).toBe('Meta.X');
  });

  it('isConfig detects config documents', () => {
    const d = new DocumentObj(makeDoc({ type: 'TestProject__config__' }), project);
    expect(d.isConfig).toBe(true);
  });

  it('isConfig is false for regular documents', () => {
    const d = new DocumentObj(makeDoc({ type: 'SomeType' }), project);
    expect(d.isConfig).toBe(false);
  });

  it('toolkit returns desc.toolkit', () => {
    const d = new DocumentObj(makeDoc({ desc: { toolkit: 'LSM' } }), project);
    expect(d.toolkit).toBe('LSM');
  });

  it('extDesc merges desc with type', () => {
    const d = new DocumentObj(makeDoc({ desc: { datasourceName: 'x' }, type: 'MyType' }), project);
    expect(d.extDesc).toEqual({ datasourceName: 'x', type: 'MyType' });
  });
});
