import numpy as np
from management import variable
from numerics_fundamentals.discretization import kgrid
from input.enum_bdry import BdryType

class user_data(object):
    def __init__(self):
        self.nspec = 8

ud = user_data()
new_state = variable.Var(5,ud)

inx = 5
iny = 1
inz = 1
x0 = 0.
x1 = 10.
y0 = 0.
y1 = 10.
z0 = 0.
z1 = 10.
left = BdryType.PERIODIC
right = BdryType.PERIODIC
bottom = BdryType.PERIODIC
top = BdryType.PERIODIC
back = BdryType.PERIODIC
front = BdryType.PERIODIC
grid = kgrid.Grid(inx,iny,inz,x0,x1,y0,y1,z0,z1,left,right,bottom,top,back,front)

print(new_state.X)