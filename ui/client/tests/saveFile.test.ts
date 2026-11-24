import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { saveFile } from '../src/utils/saveFile';

// saveFile clicks a download link and frees the bytes only later.

let revokeObjectURL: ReturnType<typeof vi.fn>;
let clickedDownloads: string[];

beforeEach(() => {
  vi.useFakeTimers();
  URL.createObjectURL = vi.fn().mockReturnValue('blob:zip') as any;
  revokeObjectURL = vi.fn();
  URL.revokeObjectURL = revokeObjectURL as any;
  clickedDownloads = [];
  vi.spyOn(HTMLAnchorElement.prototype, 'click').mockImplementation(function (this: HTMLAnchorElement) {
    clickedDownloads.push(this.download);
  });
});

afterEach(() => {
  vi.useRealTimers();
  vi.restoreAllMocks();
});

describe('saveFile', () => {
  it('clicks a link that downloads the bytes under the file name', () => {
    saveFile(new Uint8Array([1, 2, 3]).buffer, 'projects.zip');

    expect(clickedDownloads).toEqual(['projects.zip']);
  });

  it('frees the bytes only after a delay', () => {
    saveFile(new Uint8Array([1]).buffer, 'projects.zip');

    expect(revokeObjectURL).not.toHaveBeenCalled();
    vi.advanceTimersByTime(10000);
    expect(revokeObjectURL).toHaveBeenCalledWith('blob:zip');
  });
});
