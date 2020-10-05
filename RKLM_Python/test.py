# unit tester for DA infrastructure
import os
import sys, getopt
import argparse

def main ():
    parser = argparse.ArgumentParser(description='Test script')
    parser.add_argument('-N',action='store',dest='N',help='<Optional> Set ensemble size, if none is given N=1 is used.',required=False,type=int)
    parser.add_argument('-ic','--initial_conditions',action='store',dest='ic',help='<Required> Set initial conditions',required=True,choices={'aw','tv','tv_3d','rb','igw'})
    args = parser.parse_args() # collect cmd line args
    ic = args.ic

    if ic == 'bi':
        from inputs.baroclinic_instability_periodic import UserData, sol_init
    elif ic == 'tv' or ic == 'tv_2d':
        from inputs.travelling_vortex_2D import UserData, sol_init

    if args.N is None:
        N = 1
    else:
        N = args.N
    print(sol_init)
    return N, ic

if __name__ == '__main__':
    N,ic = main()
    print(N, ic)

