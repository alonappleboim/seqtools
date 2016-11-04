import os
import sys

project_dir = os.path.sep.join(sys.modules[__name__].__file__.split(os.path.sep)[:-2])
sys.path.append(project_dir)

from common.config import *

RNA_DATA_PATH = DATA_PATH + os.path.sep + 'RNA'