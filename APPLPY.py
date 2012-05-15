import sys,os
sympy_path='/home/matt/sympy'
sys.path.append(sympy_path)
os.chdir(sympy_path)
from sympy import *
x=Symbol('x');y=Symbol('y');z=Symbol('z');t=Symbol('t')
sys.displayhook=pprint
from applpy.applpy import *
print 'WELCOME TO APPLPY'
print '_________________'
print ""
Menu()
