/// <reference types="node" />
import * as fs from 'fs';
import { dirname } from 'path';
import { describe, it, expect, beforeAll, afterAll, beforeEach, afterEach, vi } from 'vitest';
import { screen, waitFor, cleanup, act, fireEvent } from '@testing-library/react';
import { createProjectViaUI, renderApp, resetStore } from './integHelpers';
import { SHARED_SERVER_URL } from './mockFactories';

vi.mock('../../src/shared/baseurl', async () => (await import('./mockFactories')).createBaseurlMock());
vi.mock('../../src/stores/useServerConstants', async () => (await import('./mockFactories')).createServerConstantsMock());
vi.mock('../../src/io/snackbar', async () => (await import('./mockFactories')).createSnackbarMock());

const REPO_DIR = '/tmp/hera-test-central-repos';

// The server runs on this machine (see globalSetup.ts), so it sees these files.
const writeFileOnServer = async (path: string, content: string) => {
  fs.mkdirSync(dirname(path), { recursive: true });
  fs.writeFileSync(path, content);
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

// The route the central folder uses: the repository shape rule lives behind it.
const jsonsInFolder = async (folder: string, withContents = false) => {
  const r = await fetch(
    `${SHARED_SERVER_URL}/repository/jsonsInFolder?folder=${encodeURIComponent(folder)}&withContents=${withContents}`);
  return JSON.parse(await r.text());
};

const getRepository = async (path: string) => {
  const r = await fetch(`${SHARED_SERVER_URL}/repository/get?path=${encodeURIComponent(path)}`);
  return { status: r.status, body: JSON.parse(await r.text()) };
};

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
    fs.rmSync(REPO_DIR, { recursive: true, force: true });
    cleanup();
  }, 15000);

  it('the folder route accepts repository files and rejects the rest', async () => {
    const { files } = await jsonsInFolder(REPO_DIR);

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

  it('an empty JSON object is not a repository', async () => {
    await writeFileOnServer(`${REPO_DIR}/empty.json`, '{}');

    const { files } = await jsonsInFolder(REPO_DIR);

    expect(files).not.toContain(`${REPO_DIR}/empty.json`);
  }, 15000);

  it('the section names may be in any case', async () => {
    await writeFileOnServer(`${REPO_DIR}/mixcase.json`, JSON.stringify({
      MyToolkit: { Config: {}, Measurements: {} },
    }));

    const { files } = await jsonsInFolder(REPO_DIR);

    expect(files).toContain(`${REPO_DIR}/mixcase.json`);
  }, 15000);

  it('the contents come with the files when asked', async () => {
    const { files, repoJsons } = await jsonsInFolder(REPO_DIR, true);

    expect(Object.keys(repoJsons).sort()).toEqual([...files].sort());
    expect(repoJsons[`${REPO_DIR}/valid.json`]).toEqual(JSON.parse(validRepo));
  }, 15000);

  it('a folder that is not there is empty, not an error', async () => {
    const { files, repoJsons } = await jsonsInFolder(`${REPO_DIR}/nowhere`, true);

    expect(files).toEqual([]);
    expect(repoJsons).toEqual({});
  }, 15000);

  it('one repository file comes with its paths resolved', async () => {
    const answer = await getRepository(`${REPO_DIR}/valid.json`);

    expect(answer.status).toBe(200);
    expect(answer.body.json).toEqual(JSON.parse(validRepo));
    expect(answer.body.resolved).toEqual(JSON.parse(validRepo));
  }, 15000);

  it('a file that is not a repository is not found', async () => {
    const answer = await getRepository(`${REPO_DIR}/broken.json`);

    expect(answer.status).toBe(404);
  }, 15000);

  it('central repo folder renders filtered files via the UI', async () => {
    localStorage.setItem('hera-central-repo-folder', REPO_DIR);
    renderApp('/CentralTestProject');

    // "Repositories" is a plain section heading now, not a tree node, so the
    // central repo folder is visible without expanding anything first.
    const folderLabel = await screen.findByText(REPO_DIR, {}, { timeout: 15000 });
    const folderItem = folderLabel.closest('.MuiTreeItem-root')!;
    const expandIcon = folderItem.querySelector('.MuiTreeItem-iconContainer')!;
    await act(async () => { fireEvent.click(expandIcon); });

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
