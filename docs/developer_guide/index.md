# Developer Guide

This section is for developers working on or extending the Hera platform. It covers architecture, data model internals, testing, and how to build new toolkits.

---

## Architecture

- **[Core Concepts](../architecture/core_concepts.md)** - Deep dive into `Project`, `ToolkitHome`, and `abstractToolkit` with class diagrams
- **[Data Layer](../architecture/data_layer.md)** - MongoDB document model, `datatypes`, and the repository JSON pipeline

## Toolkit Implementation

- **[Measurements](measurements.md)** - GIS, meteorology, and experiment toolkit internals, class hierarchies, coordinate utilities
- **[Simulations](simulations.md)** - OpenFOAM, LSM, Gaussian, wind profile, workflow system internals
- **[Risk Assessment](risk_assessment.md)** - Agent/effects/injury system, protection policies, risk area algorithms

## RAG Search

- **[RAG Search Toolkit](rag.md)** - Architecture, indexer, search API, REST server, file watcher, MkDocs plugin, CLI reference

## API Reference

- **[API Reference](api/index.md)** - Auto-generated documentation from Python source code covering all public classes and functions

## Testing

- **[Getting Started with Tests](../testing/getting_started.md)** - Setting up the test environment and running tests
- **[Test Flow](../testing/flow.md)** - Pytest session lifecycle, fixtures, comparison helpers, and the expected-output mechanism

## Roadmap

- **[Roadmap](roadmap.md)** - Planned improvements: contract-first typed interfaces, unified toolkit registry, env var DB config

## Reference

- **[Repository Schema](../reference/repository_schema.md)** - Detailed schema documentation for the data repository
- **[Environment Variables](../configuration/env_vars.md)** - Complete reference of all environment variables
