project = u"gempipe"
author = u"Gioele Lazzari"
copyright = u"Gioele Lazzari and VituloLab"


# sphinx settings: 
extensions = [
    "myst_nb",
    'sphinx_rtd_theme',
    "autoapi.extension",
    "sphinx.ext.napoleon",  # enable the google-style docstrings
    "sphinx.ext.viewcode",  # add links to highlighted source code
]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_static_path = ['static']
html_css_files = ['custom.css',]  # these need to be inside 'html_static_path'
html_logo = "static/icon_proposal.png"
html_theme_options = {'logo_only': True,}


# notebook settings: 
html_theme = "sphinx_rtd_theme"
myst_enable_extensions = [
    "dollarmath",  # parsing of dollar $ and $$ encapsulated math
    "amsmath"  # direct parsing of amsmath LaTeX environments
]
nb_execution_mode = "off"  # "force"
myst_dmath_allow_space = True  #  inline math parsed also if there are initial/final spaces 
nb_execution_timeout = 60*15  # minutes max execution seconds before timeout exception.


# autoapi settings
autoapi_add_toctree_entry = False
autoapi_dirs = ["../src/gempipe/interface"]
autoapi_ignore = ['*-checkpoint.py', '*gempipe/assets/*', '*gempipe/interface/*_utils.py']





