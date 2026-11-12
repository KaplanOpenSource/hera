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

// What a successful answer carries, and so how to read it.
enum AnswerKind {
  Json = 'application/json',
  Bytes = 'application/octet-stream',
}

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

  // One call to a named server entrypoint. Returns its answer, or null when
  // the call failed. Errors always come back as JSON.
  private async call(path: string, init: RequestInit, options: CallOptions, kind = AnswerKind.Json): Promise<any | null> {
    let ticket: Ticket | null = null;
    if (!options.silent && this.reporter) {
      ticket = this.reporter.started(options.label);
    }
    try {
      // no-store: the answer must never come from the browser cache.
      const response = await fetch(`${BASEURL}${path}`, { cache: 'no-store', ...init });
      const contentType = response.headers.get('content-type');
      if (!response.ok && contentType?.includes(AnswerKind.Json)) {
        const body = await response.json();
        this.report(options, String(body?.detail ?? body?.error ?? response.status));
        return null;
      }
      // A server without this route answers with the SPA's HTML, not an error.
      if (!response.ok || !contentType?.includes(kind)) {
        this.report(options, `${response.status} ${contentType}`);
        return null;
      }
      if (kind === AnswerKind.Bytes) {
        return await response.arrayBuffer();
      }
      return await response.json();
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

  private withQuery(path: string, params: QueryParams): string {
    const query = new URLSearchParams(params).toString();
    let suffix = '';
    if (query) {
      suffix = `?${query}`;
    }
    return `${path}${suffix}`;
  }

  protected async getJson(path: string, params: QueryParams, options: CallOptions): Promise<any | null> {
    return await this.call(this.withQuery(path, params), {}, options);
  }

  // A GET whose answer is raw bytes, not JSON.
  protected async getBytes(path: string, params: QueryParams, options: CallOptions): Promise<ArrayBuffer | null> {
    return await this.call(this.withQuery(path, params), {}, options, AnswerKind.Bytes);
  }

  protected async postJson(path: string, body: object, options: CallOptions): Promise<any | null> {
    return await this.call(path, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body),
    }, options);
  }
}
