#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 16:58:13 2020

@author: chen
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 15:55:59 2020

@author: chen
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 13:48:18 2020

@author: chen
"""


from dolfin import *
import numpy as np
from math import *
import os, sys
import matplotlib.pyplot as plt
import sympy as sp


def DS(N,r):
    # define symbol expressions of u,p,f using sympy
    x,y = sp.symbols("x[0],x[1]")
    u1 = (sp.cos(2*sp.pi*x)-1)*(sp.cos(2*sp.pi*y)-1)
    u2 =  sp.sin(2*sp.pi*x)*sp.sin(2*sp.pi*y)
    sigma = sp.diff(u2,x,1)-sp.diff(u1,y,1)
    f1 = -sp.diff(u1,x,2) - sp.diff(u1,y,2)
    f2 = -sp.diff(u2,x,2) - sp.diff(u2,y,2)
    # transform the expressions to ccode style
    u1 = sp.printing.ccode(u1)
    u2 = sp.printing.ccode(u2)
    sigma = sp.printing.ccode(sigma)
    f1 = sp.printing.ccode(f1)
    f2 = sp.printing.ccode(f2)

    # Load mesh
    mesh = RectangleMesh(Point(-1,-1), Point(1,1), N, N)
    
    # Build function space Taylor-Hood2
    V = FiniteElement("BDM", mesh.ufl_cell(), r-1)
    Q = FiniteElement("Lagrange", mesh.ufl_cell(), r) 
           
    W = FunctionSpace(mesh, V*Q)
      
    # Boundary condition
    def Dirichlt_boundary(x,on_boundary):
        return x[0] < -1+DOLFIN_EPS or x[0] > 1-DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS or x[1] < -1+DOLFIN_EPS
    # g_D
    g_D = Expression(("0","0"),degree=4)
    
    bcD = DirichletBC(W.sub(0), g_D, Dirichlt_boundary)

    
    f = Expression((f1,f2), degree=6)
    
    # Construct function space
    (u, p) = TrialFunctions(W)
    (v, q) = TestFunctions(W)
    

    a = p*q*dx - dot(curl(q),u)*dx + dot(curl(p),v)*dx + div(u)*div(v)*dx
    L = dot(f, v)*dx 
    
    # Solve
    U = Function(W)
    solve(a == L, U,bcD)
      
    # Get sub-functions
    u, p = U.split()

    # Error
    exactu = Expression((u1,u2), degree=6)
    exactp = Expression(sigma, degree=6)
    
    erru_L2 = errornorm(exactu, u, "L2")
    erru_Hdiv0 = errornorm(exactu, u, "Hdiv0")
    
    errp_L2 = errornorm(exactp, p, "L2")
    erru_Hcurl0 = errornorm(exactp, p, "Hcurl0")


    return erru_L2,erru_Hdiv0,errp_L2,erru_Hcurl0 


# Iterations and computing errors
N = np.array([2,4, 8, 16, 32, 64])
h = 1/N


erru_L2   = np.zeros((N.size, 1))
rateu_L2  = np.zeros((N.size, 1))

erru_Hdiv0   = np.zeros((N.size, 1))
rateu_Hdiv0  = np.zeros((N.size, 1))

errp_L2   = np.zeros((N.size, 1))
ratep_L2  = np.zeros((N.size, 1))

errp_Hcurl0   = np.zeros((N.size, 1))
ratep_Hcurl0  = np.zeros((N.size, 1))

r = 4
for i in range(N.size):
          errs = DS(N[i],r)
          erru_L2[i] = errs[0]
          erru_Hdiv0[i] = errs[1]
          errp_L2[i] = errs[2]
          errp_Hcurl0[i] = errs[3]
          if (i!=0):
               rateu_L2[i] = -log(erru_L2[i]/erru_L2[i-1])/log(2)
               rateu_Hdiv0[i] = -log(erru_Hdiv0[i]/erru_Hdiv0[i-1])/log(2)
               ratep_L2[i] = -log(errp_L2[i]/errp_L2[i-1])/log(2)
               ratep_Hcurl0[i] = -log(errp_Hcurl0[i]/errp_Hcurl0[i-1])/log(2)
               



print("\t\tThe error of Lagrange-BDM r=",r)       
print("h\t\t error-u-L2  rate\t error-u-Hdiv0  rate\t error-p-L2  rate  \t error-p-Hcurl0  rate")
i = 0
print("%.2e"%h[i],"&", "\t%.2e"%erru_L2[i],"&", " ---&", "\t%.2e"%erru_Hdiv0[i],"&", " ---&", "\t%.2e"%errp_L2[i],"&"," ---&","\t%.2e"%errp_Hcurl0[i],"& ---\\\\")
for i in range(1,N.size):
    print("%.2e"%h[i],"&", "\t%.2e"%erru_L2[i],"&", " %.2f"%rateu_L2[i],"&", "\t%.2e"%erru_Hdiv0[i],"&", "%.2f"%rateu_Hdiv0[i],"&", "\t%.2e"%errp_L2[i],"&","%.2f"%ratep_L2[i],"&","\t%.2e"%errp_Hcurl0[i],"&","%.2f"%ratep_Hcurl0[i],"\\\\")
print("\n")
              