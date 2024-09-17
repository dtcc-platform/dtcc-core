from . import dtcc_pb2 as proto

# Import submodules
from . import object
from . import geometry
from . import values

# Import all classes from submodules
from .object import *
from .geometry import *
from .values import *


# Collect __all__ from submodules
modules = [object, geometry, values]
__all__ = []
for module in modules:
    for name in module.__all__:
        globals()[name] = getattr(module, name)
    __all__ += module.__all__

# Additional global imports
__all__ += ["proto"]
