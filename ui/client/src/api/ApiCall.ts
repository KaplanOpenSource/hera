import { BASEURL } from '../shared/baseurl';
import { WARMING_UP } from '../shared/constants';
import { Reporter, Ticket } from './Reporter';

const SHORT_ERROR_MAX = 120;

export type CallOptions = {
  // Shown in the running toast, and in front of the error if the call fails.
  label: string,
  // No toasts at all, for background calls like the periodic reload.
  silent?: boolean,
};

export type QueryParams = { [key: string]: string };

const shorten = (text: string) => {
  if (text.length <= SHORT_ERROR_MAX) {
    return text;
  }
  return text.slice(0, SHORT_ERROR_MAX - 1) + '…';
};

// The base of every api class: the one place in the client that fetches.
export class ApiCall {
  constructor(protected reporter?: Reporter) { }

  setReporter(reporter: Reporter) {
    this.reporter = reporter;
  }

  private report(options: CallOptions, problem: string) {
    console.error(`${options.label}: ${problem}`);
    // Stay quiet while hera warms up; the startup gate handles that wait.
    if (!options.silent && problem !== WARMING_UP) {
      this.reporter?.failed(`${options.label}: ${shorten(problem)}`);
    }
  }

  // One call to a named server entrypoint. Returns its JSON answer, or null
  // when the call failed.
  private async call(path: string, init: RequestInit, options: CallOptions): Promise<any | null> {
    let ticket: Ticket | null = null;
    if (!options.silent && this.reporter) {
      ticket = this.reporter.started(options.label);
    }
    try {
      // no-store: the answer must never come from the browser cache.
      const response = await fetch(`${BASEURL}${path}`, { cache: 'no-store', ...init });
      const contentType = response.headers.get('content-type');
      // A server without this route answers with the SPA's HTML, not an error.
      if (!contentType?.includes('application/json')) {
        this.report(options, `${response.status} ${contentType}`);
        return null;
      }
      const body = await response.json();
      if (!response.ok) {
        this.report(options, String(body?.detail ?? body?.error ?? response.status));
        return null;
      }
      return body;
    } catch (e: any) {
      // The server restarting must not break the UI, so a dead connection is
      // reported like any other failure instead of throwing at the caller.
      this.report(options, `network error: ${e?.message ?? e}`);
      return null;
    } finally {
      if (ticket !== null) {
        this.reporter?.finished(ticket);
      }
    }
  }

  protected async getJson(path: string, params: QueryParams, options: CallOptions): Promise<any | null> {
    const query = new URLSearchParams(params).toString();
    let suffix = '';
    if (query) {
      suffix = `?${query}`;
    }
    return await this.call(`${path}${suffix}`, {}, options);
  }

  protected async postJson(path: string, body: object, options: CallOptions): Promise<any | null> {
    return await this.call(path, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body),
    }, options);
  }
}
