project = u"gempipe"
author = u"Gioele Lazzari"
copyright = u"Gioele Lazzari"
extensions = [
    "myst_nb",
    "autoapi.extension",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode", 
    "sphinxcontrib.mermaid"
]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "sphinx_rtd_theme"
autoapi_dirs = ["../src/gempipe"]
autoapi_ignore = ['*-checkpoint.py', '../scr/gempipe/assets/*']
nb_execution_mode = "force"



