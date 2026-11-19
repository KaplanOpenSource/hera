import { describe, it, expect } from 'vitest';
import { readFileSync } from 'fs';
import { dirname, join } from 'path';
import { fileURLToPath } from 'url';
import { unzipSync, zipSync } from 'fflate';
import { readProjectZip } from '../src/utils/readProjectZip';

const fixture = (name: string) => {
  return new Uint8Array(readFileSync(join(dirname(fileURLToPath(import.meta.url)), 'fixtures', name)));
};

describe('readProjectZip', () => {
  it('reads the name and folder from a zip exported without the progress bar', () => {
    expect(readProjectZip(fixture('ZipSample.zip'))).toEqual({
      name: 'ZipSample',
      filesDirectory: '/projects/ZipSample',
    });
  });

  it('reads the name and folder from a zip exported with the progress bar', () => {
    expect(readProjectZip(fixture('ZipSampleBatched.zip'))).toEqual({
      name: 'ZipSample',
      filesDirectory: '/projects/ZipSample',
    });
  });

  it('stops at the entry with the settings record and leaves later entries unread', () => {
    const sample = readProjectZip(fixture('ZipSample.zip'));
    const withBrokenLater = zipSync({
      ...unzipSync(fixture('ZipSample.zip')),
      'chunk_9': new TextEncoder().encode('not a pickle'),
    });
    expect(readProjectZip(withBrokenLater)).toEqual(sample);
  });

  it('reads the settings record from a later entry', () => {
    const entries = unzipSync(fixture('ZipSample.zip'));
    const settingsLast = zipSync({ 'chunk_1': entries['chunk_1'], 'chunk_0': entries['chunk_0'] });
    expect(readProjectZip(settingsLast)).toEqual({
      name: 'ZipSample',
      filesDirectory: '/projects/ZipSample',
    });
  });

  it('gives the name but no folder when the zip has no settings record', () => {
    const noSettings = zipSync({ 'chunk_1': unzipSync(fixture('ZipSample.zip'))['chunk_1'] });
    expect(readProjectZip(noSettings)).toEqual({ name: 'ZipSample', filesDirectory: '' });
  });

  it('throws on a zip that holds no pickled records', () => {
    const notAnExport = zipSync({ 'readme.txt': new TextEncoder().encode('hello') });
    expect(() => readProjectZip(notAnExport)).toThrow();
  });
});
