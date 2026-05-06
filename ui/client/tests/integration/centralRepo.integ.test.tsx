/// <reference types="node" />
import { describe, it, expect, beforeAll, afterAll, beforeEach, afterEach, vi } from 'vitest';
import { screen, waitFor, cleanup, act, fireEvent } from '@testing-library/react';
import { createProjectViaUI, renderApp, resetStore } from './integHelpers';
import { SHARED_SERVER_URL } from './mockFactories';

vi.mock('../../src/shared/baseurl', async () => (await import('./mockFactories')).createBaseurlMock());
vi.mock('../../src/stores/useServerConstants', async () => (await import('./mockFactories')).createServerConstantsMock());
vi.mock('../../src/io/snackbar', async () => (await import('./mockFactories')).createSnackbarMock());

const execPython = async (code: string) => {
  const r = await fetch(`${SHARED_SERVER_URL}/exec`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ code }),
  });
  return JSON.parse(await r.text());
};

const REPO_DIR = '/tmp/hera-test-central-repos';

const writeFileOnServer = async (path: string, content: string) => {
  const escaped = content.replace(/\\/g, '\\\\').replace(/'/g, "\\'").replace(/\n/g, '\\n');
  await execPython(`
import os
path = '${path}'
os.makedirs(os.path.dirname(path), exist_ok=True)
with open(path, 'w') as f:
    f.write('${escaped}')
result = 'ok'
`);
};

const validRepo = JSON.stringify({
  MyToolkit: {
    config: {},
    measurements: {},
  },
});

const validRepoMultiToolkit = JSON.stringify({
  ToolkitA: { config: {}, datasource: {} },
  ToolkitB: { simulations: {}, cache: {} },
});

const invalidNotDict = JSON.stringify([1, 2, 3]);

const invalidBadSection = JSON.stringify({
  MyToolkit: {
    config: {},
    unknownSection: {},
  },
});

const invalidToolkitNotDict = JSON.stringify({
  MyToolkit: "not a dict",
});

const notJson = 'this is not json {{{';

describe('Central Repo isRepoJson integration', () => {
  beforeAll(async () => {
    await createProjectViaUI('CentralTestProject');

    await writeFileOnServer(`${REPO_DIR}/valid.json`, validRepo);
    await writeFileOnServer(`${REPO_DIR}/multi.json`, validRepoMultiToolkit);
    await writeFileOnServer(`${REPO_DIR}/sub/nested.json`, validRepo);
    await writeFileOnServer(`${REPO_DIR}/notdict.json`, invalidNotDict);
    await writeFileOnServer(`${REPO_DIR}/badsection.json`, invalidBadSection);
    await writeFileOnServer(`${REPO_DIR}/toolkitnotdict.json`, invalidToolkitNotDict);
    await writeFileOnServer(`${REPO_DIR}/broken.json`, notJson);
    await writeFileOnServer(`${REPO_DIR}/caseConfiguration.json`, validRepo);
    await writeFileOnServer(`${REPO_DIR}/readme.txt`, 'not a json file');
  }, 90000);

  beforeEach(() => { resetStore(); });
  afterEach(() => { cleanup(); });

  afterAll(async () => {
    await execPython(`
import shutil
shutil.rmtree('${REPO_DIR}', ignore_errors=True)
result = 'ok'
`);
    cleanup();
  }, 15000);

  it('isRepoJson accepts valid repo and rejects invalid files', async () => {
    const response = await execPython(`
import os, glob, json

SECTIONS = {'config', 'datasource', 'measurements', 'simulations', 'cache', 'function'}

def isRepoJson(path):
    try:
        with open(path, 'r') as f:
            doc = json.load(f)
    except Exception:
        return False
    if not isinstance(doc, dict) or not doc:
        return False
    for toolkitValue in doc.values():
        if not isinstance(toolkitValue, dict):
            return False
        for sectionKey in toolkitValue.keys():
            if sectionKey.lower() not in SECTIONS:
                return False
    return True

folder = '${REPO_DIR}'
allFiles = glob.glob(os.path.join(folder, '**', '*.json'), recursive=True)
allFiles = [f for f in allFiles if not f.endswith('caseConfiguration.json')]
jsonFiles = sorted(f for f in allFiles if isRepoJson(f))
result = jsonFiles
`);

    const files: string[] = response.data;

    expect(files).toContain(`${REPO_DIR}/multi.json`);
    expect(files).toContain(`${REPO_DIR}/valid.json`);
    expect(files).toContain(`${REPO_DIR}/sub/nested.json`);

    expect(files).not.toContain(`${REPO_DIR}/notdict.json`);
    expect(files).not.toContain(`${REPO_DIR}/badsection.json`);
    expect(files).not.toContain(`${REPO_DIR}/toolkitnotdict.json`);
    expect(files).not.toContain(`${REPO_DIR}/broken.json`);
    expect(files).not.toContain(`${REPO_DIR}/caseConfiguration.json`);

    expect(files).toHaveLength(3);
  }, 30000);

  it('isRepoJson rejects empty JSON object', async () => {
    await writeFileOnServer(`${REPO_DIR}/empty.json`, '{}');

    const response = await execPython(`
import json

SECTIONS = {'config', 'datasource', 'measurements', 'simulations', 'cache', 'function'}

def isRepoJson(path):
    try:
        with open(path, 'r') as f:
            doc = json.load(f)
    except Exception:
        return False
    if not isinstance(doc, dict) or not doc:
        return False
    for toolkitValue in doc.values():
        if not isinstance(toolkitValue, dict):
            return False
        for sectionKey in toolkitValue.keys():
            if sectionKey.lower() not in SECTIONS:
                return False
    return True

result = isRepoJson('${REPO_DIR}/empty.json')
`);

    expect(response.data).toBe(false);
  }, 15000);

  it('isRepoJson is case-insensitive on section keys', async () => {
    await writeFileOnServer(`${REPO_DIR}/mixcase.json`, JSON.stringify({
      MyToolkit: { Config: {}, Measurements: {} },
    }));

    const response = await execPython(`
import json

SECTIONS = {'config', 'datasource', 'measurements', 'simulations', 'cache', 'function'}

def isRepoJson(path):
    try:
        with open(path, 'r') as f:
            doc = json.load(f)
    except Exception:
        return False
    if not isinstance(doc, dict) or not doc:
        return False
    for toolkitValue in doc.values():
        if not isinstance(toolkitValue, dict):
            return False
        for sectionKey in toolkitValue.keys():
            if sectionKey.lower() not in SECTIONS:
                return False
    return True

result = isRepoJson('${REPO_DIR}/mixcase.json')
`);

    expect(response.data).toBe(true);
  }, 15000);

  it('central repo folder renders filtered files via the UI', async () => {
    localStorage.setItem('hera-central-repo-folder', REPO_DIR);
    renderApp('/CentralTestProject');

    const reposLabel = await screen.findByText('Repositories', {}, { timeout: 15000 });
    const reposContent = reposLabel.closest('.MuiTreeItem-content')!;
    await act(async () => { fireEvent.click(reposContent); });

    const folderLabel = await screen.findByText(REPO_DIR, {}, { timeout: 5000 });
    const centralContent = folderLabel.closest('.MuiTreeItem-content')!;
    await act(async () => { fireEvent.click(centralContent); });

    await waitFor(() => {
      expect(screen.getByText(`${REPO_DIR}/valid.json`)).toBeTruthy();
      expect(screen.getByText(`${REPO_DIR}/multi.json`)).toBeTruthy();
      expect(screen.getByText(`${REPO_DIR}/sub/nested.json`)).toBeTruthy();
    }, { timeout: 15000 });

    expect(screen.queryByText(`${REPO_DIR}/broken.json`)).toBeNull();
    expect(screen.queryByText(`${REPO_DIR}/badsection.json`)).toBeNull();
    expect(screen.queryByText(`${REPO_DIR}/notdict.json`)).toBeNull();
    expect(screen.queryByText(`${REPO_DIR}/caseConfiguration.json`)).toBeNull();
  }, 30000);
}, 180000);
