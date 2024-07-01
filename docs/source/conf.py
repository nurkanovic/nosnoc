# Configuration file for the Sphinx documentation builder.

import os

# -- Project information

project = 'nosnoc'
copyright = '2024, Armin Nurkanovic, Jonathan Frey, Anton Pozharskiy, Moritz Diehl'
author = 'Armin Nurkanovic'

release = '0.5'
version = '0.5.0'

# -- General configuration


extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinxcontrib.matlab',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx_math_dollar',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

matlab_src_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
matlab_keep_package_prefix = False
primary_domain = "mat"
default_role = "obj"
autodoc_member_order = 'bysource'

napoleon_numpy_docstring = False
napoleon_preprocess_types = True
napoleon_use_param = False
napoleon_attr_annotations = False

autodoc_default_options = {'members': True, 'show-inheritance': True}
autosummary_generate = True

# -- magic to kill mathjax

mathjax3_config = {
  "tex": {
    "inlineMath": [['\\(', '\\)']],
    "displayMath": [["\\[", "\\]"]],
  }
}
