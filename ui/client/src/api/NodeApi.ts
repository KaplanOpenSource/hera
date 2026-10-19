import { ApiCall } from './ApiCall';
import { NodeCatalogEntry } from '../components/workflow/nodeCatalog';

export interface NodeParamsValidation {
  ok: boolean;
  // Hermes's own message for the first bad parameter. Empty when there is
  // nothing to say (all good, unknown type, or a node with no checks).
  message: string;
}

// Every call to /node-catalog and /node-params. Each method answers data, or null when the call
// failed.
export class NodeApi extends ApiCall {
  // Tests one node's parameter values with the node's own Hermes checks. The
  // server reaches hera and MongoDB for this, so call it sparingly (debounced).
  async validateParams({
    type,
    params,
  }: {
    type: string,
    params: { [param: string]: unknown },
  }): Promise<NodeParamsValidation | null> {
    return await this.postJson('/node-params/validate', { type, params }, { label: 'node params validation' });
  }

  // The Hermes node catalog: available node types and their parameters/outputs.
  async catalog(): Promise<NodeCatalogEntry[] | null> {
    return await this.getJson('/node-catalog', {}, { label: 'node catalog' });
  }
}
