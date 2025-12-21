"""
gmpy2 documentation build configuration file.

This file is execfile()d with the current directory set to its containing dir.

The contents of this file are pickled, so don't put values in the namespace
that aren't pickleable (module imports are okay, they're removed
automatically).
"""

import os
import sys

import packaging.version

import gmpy2

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.imgmath', 'sphinx.ext.doctest',
              'sphinx.ext.intersphinx', 'sphinx_rtd_theme']

# The name of a reST role (builtin or Sphinx extension) to use as the
# default role, that is, for text marked up `like this`.
default_role = 'py:obj'

# Sphinx will warn about all references where the target cannot be found.
nitpicky = True

# This value selects if automatically documented members are sorted
# alphabetical (value 'alphabetical'), by member type (value 'groupwise')
# or by source order (value 'bysource').
autodoc_member_order = 'groupwise'

# The default options for autodoc directives. They are applied to all
# autodoc directives automatically.
autodoc_default_options = {'members': True}

# Contains mapping the locations and names of other projects that
# should be linked to in this documentation.
intersphinx_mapping = {'python': ('https://docs.python.org/3/', None)}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# General information about the project.
project = gmpy2.__package__
copyright = '2012 - 2025, Case Van Horsen'

gmpy2_version = packaging.version.parse(gmpy2.__version__)

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = f"{gmpy2_version.major}.{gmpy2_version.minor}"
# The full version, including alpha/beta/rc tags.
release = gmpy2.__version__

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'sphinx_rtd_theme'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [('index', 'gmpy2.tex', 'gmpy2 Documentation',
                    'Case Van Horsen', 'manual')]

# A dictionary that contains LaTeX snippets overriding those Sphinx
# usually puts into the generated .tex files.
latex_elements = {
    'preamble': r'''
\DeclareUnicodeCharacter{03B5}{$\epsilon$}
''',
}

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [('index', 'gmpy2', 'gmpy2 Documentation', ['Case Van Horsen'], 3)]

# Python code that is treated like it were put in a testcleanup directive
# for *every* file that is tested, and for every group.
doctest_global_cleanup = """
import gmpy2

gmpy2.set_context(gmpy2.context())
"""
