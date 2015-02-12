# -*- coding: utf-8 -*-
"""
fmincon alternative using CasADi and IpOpt:

Implement your own objective funtion, gradient and hessian. 
You need to use CasADi's data structures, and here is guide: ../Resources/casadi-users_guide.pdf
Then use fmincon(N, x0, lb, ub, TolX, MaxIter) to solve it, as shown in main thread at the bottom.

A test case is included in this file.
It's Matlab version is in ../fmincon_alternative/testcase.m

To do:
1. Integrate fmincon piece into our python software. To do this:
    * install CasADi and IpOpt. Please look into Python/README
    * implement set_Q_dV_unc in Python and use fmincon() to solve it. Unlike matlab version, in Python you have to write them separately.
        ** set_Q_dV_unc
        ** gradient of set_Q_dV_unc
        ** hessian of set_Q_dV_unc
    * you need to use CasADi’s data structures in set_Q_dV_unc. Please look into Resources/casadi-users_guide.pdf
2. verify the correctness of fmincon alternative. Though Yunye and I tested this function using a small test case (this test case is included in fmincon.py, and fmincon_alternative/testcase.m), it could happen that it goes to something wrong with our software’s case.

by Yunhan Wang, yw559
Updated Sep 29, 2014
"""

import numpy as np
from casadi import *

# objective function
def objF(l):
    """ test case, implement your own objective function """
    return 1/(4*pi**2*l[0]*l[1])*exp(-0.5*((x1[0]**2+x1[1]**2)/l[0]+(x2[0]**2+x2[1]**2)/l[1]))

# gradient function
def gradF(l):
    """ test case, implement your own gradient function """
    q= 1/(4*pi**2*l[0]*l[1])*exp(-0.5*((x1[0]**2+x1[1]**2)/l[0]+(x2[0]**2+x2[1]**2)/l[1]))
    gq=MX.zeros(N)
    gq[0]=q*(x1[0]**2+x1[1]**2-2*l[0])/(2*l[0]**2)
    gq[1]=q*(x2[0]**2+x2[1]**2-2*l[1])/(2*l[1]**2)
    return gq

# hessian of language
def hessF(l):
    """ test case, implement your own hessian function """
    q= 1/(4*pi**2*l[0]*l[1])*exp(-0.5*((x1[0]**2+x1[1]**2)/l[0]+(x2[0]**2+x2[1]**2)/l[1]))
    hq=MX.zeros(N,N)
    hq[0,0]=q*((x1[0]**2+x1[1]**2)**2+8*l[0]**2-8*l[0]*(x1[0]**2+x1[1]**2))/(4*l[0]**2);
    hq[1,1]=q*((x2[0]**2+x2[1]**2)**2+8*l[1]**2-8*l[1]*(x2[0]**2+x2[1]**2))/(4*l[1]**2);
    hq[0,1]=q*(x1[0]**2+x1[1]**2-2*l[0])/(2*l[0]**2)*(x2[0]**2+x2[1]**2-2*l[1])/(2*l[1]**2)
    hq[1,0]=hq[0,1]
    return hq
    
def fmincon(N, l0, lb, ub, TolX, MaxIter):
    """
    use IpOpt (linear search & interior point) to solve NLP
    do NOT change this function if you are not familiar
    """
    # Optimization variables declaration & init guess
    l = MX.sym("x",N) # variables, sym must be called "x" in order for solver to find
    
    # consistent gradient function name & hessian function name
    grad = MXFunction(gradFIn(x=l), gradFOut(grad=gradF(l)))
    h = MXFunction(hessLagIn(x=l), hessLagOut(hess=hessF(l)))
    
    # Create NLP
    # g is s.t function that we don't need, so I give 0<=0<=0 that always holds
    nlp=MXFunction(nlpIn(x=l),nlpOut(f=objF(l),g=0)) 
    solver = NlpSolver("ipopt", nlp)

    # Options of NLP
    # use solver.printOptions() to see avaible options

    #solver.setOption("linear_solver","ma27") # HSL linear solver is said to be faster
    solver.setOption("max_iter", MaxIter)
    solver.setOption("acceptable_tol",TolX)
    #solver.setOption("print_level", 9)
    solver.setOption("hessian_approximation","exact")
    solver.setOption("grad_f", grad)
    solver.setOption("hess_lag", h)
    solver.setOption("derivative_test", "first-order")

    # Initialize solver
    solver.init()
    solver.setInput(l0, "x0")
    # Initial guess and bounds for variables
    solver.setInput(lb,"lbx") # change lower bound and upper bound here!
    solver.setInput(ub,"ubx")
    # nonlinear bounds, none. corresponding with g=0 (0<=0<=0 always holds)
    solver.setInput(0,"lbg")
    solver.setInput(0,"ubg")

    solver.evaluate()
    return solver.getOutput("x"), solver.getOutput("f")
    
if __name__ == "__main__":
    # x1, x2 are constant required by this testcase only
    x1 = DMatrix([0.8147, 0.9058])
    x2 = DMatrix([0.1270, 0.9134])
    
    N = 2 # number of variables
    x0 = [0.6324, 0.0975] # initial guess of variables
    lb = [-1]*N # lower bound of variables
    ub = [1]*N # upper bound of variables
    TolX = 1e-5 # tolerance, same as in Matlab
    MaxIter = 20 # max # of iterations to run
    
    # command to solve NLP
    newl,newq = fmincon(N, x0, lb, ub, TolX, MaxIter)
    print "nu*=",newl
    print "optimal=", newq