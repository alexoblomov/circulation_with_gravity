#driver script for paper plots
from parameters import *
from varyG_compliances import *
from varyG_heights import *
from varyG_Pthorax import *
from varyG_Vtotal import * 



with open("scripts/varyG.py") as f:
    exec(f.read())
with open("scripts/varyG_compliances.py") as f:
    exec(f.read())
with open("scripts/varyG_heights.py") as f:
    exec(f.read())
with open("scripts/varyG_Vtotal.py") as f:
    exec(f.read())
