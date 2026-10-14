import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, cleanup, waitFor } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';

vi.mock('../src/io/entrypoint', () => ({
  getJson: vi.fn(),
}));

import { getJson } from '../src/io/entrypoint';
import { FetchProjects } from '../src/io/FetchProjects';
import { useProjectStore } from '../src/stores/useProjectStore';
import { useProjectListStore } from '../src/stores/useProjectListStore';

const renderFetch = () => {
  return render(
    <MemoryRouter initialEntries={['/Alpha']}>
      <FetchProjects urlProjectName="Alpha" />
    </MemoryRouter>
  );
};

beforeEach(() => {
  vi.mocked(getJson).mockReset();
  vi.mocked(getJson).mockResolvedValue(null);
  useProjectListStore.setState({ projectNames: [{ name: 'Alpha' }, { name: 'Beta' }] });
});

afterEach(() => {
  cleanup();
});

describe('FetchProjects after startup', () => {
  it('does not fetch the project startup already loaded', async () => {
    useProjectStore.setState({
      currProjectName: 'Alpha',
      currProject: { name: 'Alpha', documents: [] },
    });

    renderFetch();

    await waitFor(() => {
      expect(useProjectStore.getState().currProjectName).toBe('Alpha');
    });
    expect(getJson).not.toHaveBeenCalled();
  });

  it('fetches when the store holds another project', async () => {
    useProjectStore.setState({
      currProjectName: 'Alpha',
      currProject: { name: 'Beta', documents: [] },
    });

    renderFetch();

    await waitFor(() => {
      expect(getJson).toHaveBeenCalled();
    });
    const paths = vi.mocked(getJson).mock.calls.map(c => c[0]);
    expect(paths).toEqual(['/project/toolkits', '/project/documents']);
  });
});
