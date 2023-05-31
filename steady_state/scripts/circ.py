#driver script for paper plots
from parameters import *
from vary_compliances import *
from vary_heights import *
from varyG import *
from varyVtotal import * 



with open("scripts/varyG.py") as f:
    exec(f.read())
with open("scripts/vary_compliances.py") as f:
    exec(f.read())
with open("scripts/vary_heights.py") as f:
    exec(f.read())
with open("scripts/varyVtotal.py") as f:
    exec(f.read())
