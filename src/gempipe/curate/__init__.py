"""
Please note that ``gempipe.curate`` functions are all called from ``gempipe``. 
Therefore, for example, the following two import statements are equivalent::

    # equivalent imports:
    from gempipe.curate import perform_gapfill
    from gempipe import perform_gapfill 
"""


from .gaps import *
from .egcs import *