"""
Please note that ``gempipe.curate`` functions are all called from ``gempipe``. 
Therefore, for example, the following two import statements are equivalent::

    # equivalent imports:
    from gempipe.curate import perform_gapfilling
    from gempipe import perform_gapfilling
"""



from .gaps import *
from .egcs import *
from .medium import *


# set up the cobra solver
cobra.Configuration().solver = "glpk_exact"