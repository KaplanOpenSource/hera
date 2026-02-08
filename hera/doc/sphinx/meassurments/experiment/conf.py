import os
import sys

# Ensure the Hera repository root is importable when building docs from:
# hera/doc/sphinx/meassurments/experiment
# This lets autodoc import `hera.*` modules.
sys.path.insert(0, os.path.abspath('../../../../..'))

# -- Project information -----------------------------------------------------
# These fields affect only the generated documentation UI (titles/footer/meta),
# and do NOT affect runtime logic of Hera.
project = 'Hera Documentation'
author = 'Hera Team'
copyright = ''

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',   # Extract docstrings from Python modules/classes/functions
    'sphinx.ext.napoleon',  # Support NumPy/Google style docstrings (Parameters/Returns/etc.)
    'sphinx.ext.viewcode',  # Add links to highlighted source code pages
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Optional: makes it easier to document class __init__ docstrings and show type hints nicely
autoclass_content = 'both'                 # include class docstring + __init__ docstring
autodoc_typehints = 'description'          # show type hints in the description, not only signatures
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_param = True
napoleon_use_rtype = True

# Avoid Sphinx import failures from optional dependencies during doc build.
# If you prefer strict builds, you can remove this list.
autodoc_mock_imports = [
    'argos',            # optional in your environment
    'tb_rest_client',   # optional
    'qgl',              # optional
    'FreeCAD',          # optional
]

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
