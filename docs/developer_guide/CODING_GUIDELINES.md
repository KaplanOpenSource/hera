# Coding Guidelines

## Purpose
This document defines the **coding style** that all Hera source code must follow. Consistent style improves readability, simplifies code reviews, and enables reliable linting across the project.

## Identifier Naming
| Element | Required style | Example |
|---------|----------------|---------|
| **Variables, functions, methods** | `snake_case` | `my_variable = 5` |
| **Modules / files** | `snake_case` (all lower) | `data_processing.py` |
| **Constants** | `UPPER_SNAKE_CASE` | `PI = 3.14` |
| **Public class names / Toolkit classes** | `PascalCase` (exception) | `MyToolkit` |
| **Private/internal helpers** | `_snake_case` (leading underscore) | `_helper_func` |

---
*This guideline is version‑controlled and will evolve as the codebase grows.*