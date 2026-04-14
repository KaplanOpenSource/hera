import { expect } from 'vitest';
import { render, screen, fireEvent, waitFor, within, cleanup, act } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';
import type { ProjectDocument } from '../../src/shared/types';
import App from '../../src/App';
import { NO_PROJECT, useProjectStore } from '../../src/stores/useProjectStore';

export const resetStore = () => {
  useProjectStore.setState({
    projectNames: [],
    currProjectName: NO_PROJECT,
    currProject: null,
    toolkits: [],
  });
};

/** Render the full app at a given path. */
export const renderApp = (path = '/') => {
  render(
    <MemoryRouter initialEntries={[path]}>
      <App />
    </MemoryRouter>,
  );
};

/** Render the app, click "Add project", fill the dialog, submit, cleanup. */
export const createProjectViaUI = async (projectName: string) => {
  resetStore();
  renderApp('/');

  const addWrapper = await screen.findByLabelText('Add project');
  await act(async () => { fireEvent.click(within(addWrapper).getByRole('button')); });

  const dialog = await screen.findByRole('dialog');
  fireEvent.change(within(dialog).getByRole('textbox', { name: /project name/i }), {
    target: { value: projectName },
  });
  const loadRepos = within(dialog).getByRole('checkbox', { name: /load repositories/i });
  if ((loadRepos as HTMLInputElement).checked) fireEvent.click(loadRepos);

  await act(async () => {
    fireEvent.click(within(dialog).getByText('Add Project'));
  });

  await waitFor(() => {
    const names = useProjectStore.getState().projectNames;
    expect(names.some(p => p.name === projectName)).toBe(true);
  }, { timeout: 15000 });

  cleanup();
};

/** Render the app with a project loaded, click "Add Document", fill the dialog, submit, cleanup. Returns the new doc id. */
export const addDocumentViaUI = async (
  projectName: string,
  docName: string,
  opts: { agent?: boolean } = {},
): Promise<string> => {
  resetStore();
  renderApp('/' + encodeURIComponent(projectName));

  await waitFor(() => {
    expect(useProjectStore.getState().currProject?.name).toBe(projectName);
  }, { timeout: 15000 });

  const addWrapper = await screen.findByLabelText('Add Document');
  await act(async () => { fireEvent.click(within(addWrapper).getByRole('button')); });

  const dialog = await screen.findByRole('dialog');
  fireEvent.change(within(dialog).getByRole('textbox', { name: /^name$/i }), {
    target: { value: docName },
  });
  if (opts.agent) {
    fireEvent.mouseDown(within(dialog).getByText('Regular'));
    fireEvent.click(await screen.findByRole('option', { name: 'Agent' }));
  }

  await act(async () => {
    fireEvent.click(within(dialog).getByRole('button', { name: /add document/i }));
  });

  await waitFor(() => {
    const docs: ProjectDocument[] = useProjectStore.getState().currProject?.documents ?? [];
    expect(docs.some(d => d.desc.datasourceName === docName)).toBe(true);
  }, { timeout: 15000 });

  const docs: ProjectDocument[] = useProjectStore.getState().currProject!.documents;
  const doc = docs.find(d => d.desc.datasourceName === docName)!;

  cleanup();
  return doc._id.$oid;
};
