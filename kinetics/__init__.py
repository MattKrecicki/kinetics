#from ntpSystem.comparison import *
from ntpSystem.containers import *
from ntpSystem.data import *
from ntpSystem.debug import *
from ntpSystem.errors import *
from ntpSystem.examples import *
from ntpSystem.functions import *
from ntpSystem.containers import *
from ntpSystem.jupyter import *
from ntpSystem.tests import *
from ntpSystem.simplepke import *


import os

_ROOT = os.path.abspath(os.path.dirname(__file__))
def setDataPath(path):
    return os.path.join(_ROOT, 'data', path)