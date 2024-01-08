project = u"gempipe"
author = u"Gioele Lazzari"
copyright = u"Gioele Lazzari"
extensions = [
    "myst_nb",
    "autoapi.extension",
    "sphinx.ext.napoleon",  # enable the google-style docstrings
    "sphinx.ext.viewcode",  # add links to highlighted source code
    "sphinx.ext.mathjax",  # view latex equations
]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "sphinx_rtd_theme"
autoapi_add_toctree_entry = False
autoapi_dirs = ["../src/gempipe/curate"]
autoapi_ignore = ['*-checkpoint.py', '*gempipe/assets/*', ]
nb_execution_mode = "force"
myst_enable_extensions = [
    "dollarmath",  # parsing of dollar $ and $$ encapsulated math
    "amsmath"  # direct parsing of amsmath LaTeX environments
]
myst_dmath_allow_space = True  #  inline math parsed also if there are initial/final spaces 
nb_execution_timeout = 60*5  # minutes max execution seconds before timeout exception.



