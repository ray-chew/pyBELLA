from discretization.kgrid import Grid, ElemSpaceDiscr, NodeSpaceDiscr
from management.variable import States

import numpy as np

def data_init(ud):
    inx = ud.inx
    iny = ud.iny
    inz = ud.inz
    x0 = ud.xmin
    x1 = ud.xmax
    y0 = ud.ymin
    y1 = ud.ymax
    z0 = ud.zmin
    z1 = ud.zmax
    left = ud.bdry_type_min[0]
    right = ud.bdry_type_max[0]
    bottom = ud.bdry_type_min[1]
    top = ud.bdry_type_max[1]
    back = ud.bdry_type_min[2]
    front = ud.bdry_type_max[2]

    grid = Grid(inx,iny,inz,x0,x1,y0,y1,z0,z1,left,right,bottom,top,back,front)

    elem = ElemSpaceDiscr(grid)
    node = NodeSpaceDiscr(grid)

    return elem, node


