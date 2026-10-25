# Developer Guide

This section is for developers working on or extending the Hera platform. It covers architecture, data model internals, testing, and how to build new toolkits.

---

## Architecture

- **[Core Concepts](../architecture/core_concepts.md)** - Deep dive into `Project`, `ToolkitHome`, `abstractToolkit`, and `dataToolkit` with class diagrams
- **[Data Layer](../architecture/data_layer.md)** - MongoDB document model, `datatypes`, the repository JSON pipeline, connections, and function caching

## Toolkit Implementation

- **[Measurements](measurements/index.md)** - GIS, meteorology, and experiment toolkit internals, class hierarchies, coordinate utilities
- **[Simulations](simulations/index.md)** - OpenFOAM, LSM, Gaussian, wind profile, ML/DL, workflow system internals
- **[Risk Assessment](risk_assessment.md)** - Agent/effects/injury system, protection policies, risk area algorithms

## RAG Search

- **[RAG Search Toolkit](rag.md)** - Architecture, indexer, search API, REST server, file watcher, MkDocs plugin, CLI reference

## API Reference

- **[API Reference](api/index.md)** - Auto-generated documentation from Python source code covering all public classes and functions

## Testing

- **[Getting Started with Tests](../testing/getting_started.md)** - Setting up the test environment and running tests
- **[Test Flow](../testing/flow.md)** - Pytest test-session lifecycle, fixtures, repository data loading, and teardown
- **[Unit Tests](../testing/unit_tests.md)** - The hermetic unit-test layer under `hera/tests/unit/` (pending merge into master)

## Roadmap

- **[Roadmap](roadmap.md)** - Planned improvements: contract-first typed interfaces, unified toolkit registry, env var DB config, async DB access, data-handler plugins

## Reference

- **[Coding Guidelines](CODING_GUIDELINES.md)** — Naming conventions and style rules for the Hera codebase (snake_case, exceptions, linting).

- **[Repository Schema](../reference/repository_schema.md)** - Detailed schema documentation for the data repository
- **[Environment Variables](../configuration/env_vars.md)** - Complete reference of all environment variables
