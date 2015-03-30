from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
sys.path.append(".")

from .variables import *
from .model_part import *
from .node import *
from .gid_io import *
from . import triangle
from . import solving_strategy
from . import static_scheme
from . import builder_and_solver
from . import line2d
