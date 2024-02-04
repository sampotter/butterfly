# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'butterfly'
copyright = '2023, Samuel F. Potter'
author = 'Samuel F. Potter'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['breathe']

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_permalinks_icon = '§'
html_theme = 'insipid'
html_static_path = ['_static']


# -- Breathe Configuration ---------------------------------------------------

breathe_projects = {'butterfly': './xml'}
breathe_default_project = 'butterfly'

# see https://github.com/breathe-doc/breathe/issues/477
breathe_domain_by_extension = {"c": "c", "h": "c", "cpp": "cpp", "hpp": "cpp"}
