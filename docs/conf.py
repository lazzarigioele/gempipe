project = u"gempipe"
author = u"Gioele Lazzari"
copyright = u"Gioele Lazzari"
extensions = [
    "myst_nb",
    "autoapi.extension",
    "sphinx.ext.napoleon",  # enable the google-style docstrings
    "sphinx.ext.viewcode",  # add links to highlighted source code
]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "sphinx_rtd_theme"
autoapi_add_toctree_entry = False
autoapi_dirs = ["../src/gempipe/curate"]
autoapi_ignore = ['*-checkpoint.py', '*gempipe/assets/*', ]
nb_execution_mode = "force"




