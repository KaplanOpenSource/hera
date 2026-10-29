import { describe, it, expect } from 'vitest';
import { readdirSync, readFileSync } from 'fs';
import { join } from 'path';

const API_DIR = join(__dirname, '..', 'src', 'api');

// The api package is every call to the server and nothing else. It must not
// reach into React, the stores, the components, or the client's toasts.
const FORBIDDEN = ['react', '/stores/', '/io/', '/components/', 'notistack'];

const apiFiles = () => {
  return readdirSync(API_DIR).filter(name => name.endsWith('.ts'));
};

describe('src/api', () => {
  it('has files', () => {
    expect(apiFiles().length).toBeGreaterThan(0);
  });

  it.each(apiFiles())('%s imports nothing from React, a store, io or a component', (name) => {
    const text = readFileSync(join(API_DIR, name), 'utf8');
    const imported: string[] = [];
    for (const line of text.split('\n')) {
      const match = line.match(/from\s+['"]([^'"]+)['"]/);
      if (match) {
        imported.push(match[1]);
      }
    }
    for (const source of imported) {
      for (const forbidden of FORBIDDEN) {
        expect(source.toLowerCase().includes(forbidden), `${name} imports ${source}`).toBe(false);
      }
    }
  });
});
