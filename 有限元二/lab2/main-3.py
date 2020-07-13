#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 11:00:10 2019

@author: chen
"""



from dolfin import *
import numpy as np
from math import *


def Stokes(N,mu,gamma):
    # Load mesh
    mesh = UnitSquareMesh(N,N)
    
    # Build function space
    V = VectorElement("Lagrange", mesh.ufl_cell(), 2)
    Q = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    W = FunctionSpace(mesh, V*Q)
    
    # Next, we define the boundary conditions. ::
    
    # Boundaries
    def Dirichlt_boundary(x,on_boundary):
        return x[0] < DOLFIN_EPS or x[0] > 1-DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS
    
    # g_D
    
    g_D = Expression(("pow(x[0],2)*pow(1-x[0],2)*x[1]*(1-x[1])*(1-2*x[1])",
                      "-x[0]*(1-x[0])*(1-2*x[0])*pow(x[1],2)*pow(1-x[1],2)"),degree=10)
    
    
    #Boundaries conditions
    bcD = DirichletBC(W.sub(0), g_D, Dirichlt_boundary)
    
    
    

    
    #Source term f
    
   
      
    f0 = Expression(("-(4*pow(x[0],2)*pow(x[0]-1,2)*(x[1] - 1) + 2*pow(x[0],2)*(2*x[1] - 1)*pow(x[0]-1,2) \
                + 4*pow(x[0],2)*x[1]*pow(x[0]-1,2) + 2*x[1]*(2*x[1] - 1)*pow(x[0]-1,2)*(x[1] - 1)\
                    + 2*pow(x[0],2)*x[1]*(2*x[1] - 1)*(x[1] - 1) + 4*x[0]*x[1]*(2*x[0] - 2)*(2*x[1] - 1)*(x[1] - 1))",\
        "(4*pow(x[1],2)*pow(x[1]-1,2)*(x[0] - 1) + 2*pow(x[1],2)*(2*x[0] - 1)*pow(x[1]-1,2) \
        + 4*pow(x[1],2)*x[0]*pow(x[1]-1,2) + 2*x[0]*(2*x[0] - 1)*pow(x[1]-1,2)*(x[0] - 1) \
            + 2*pow(x[1],2)*x[0]*(2*x[0] - 1)*(x[0] - 1) + 4*x[1]*x[0]*(2*x[1] - 2)*(2*x[0] - 1)*(x[0] - 1))\
    "),degree=6)
    f1 = Expression(("10*(3*pow(x[0]-0.5,2)*pow(x[1],2)-3*pow(1-x[0],2)*pow(x[1]-0.5,3))",\
                 "10*(2*pow(x[0]-0.5,3)*x[1]+3*pow(1-x[0],3)*pow(x[1]-0.5,2))"), degree =6)
     
    # Construct function space
    (u, p) = TrialFunctions(W)
    (v, q) = TestFunctions(W)
    
    # Define stress and strain
    def epsilon(u):
        return 0.5*(grad(u)+grad(u).T)
    
    # Define variational problem
    a = 2*mu*inner(epsilon(u), epsilon(v))*dx - div(v)*p*dx + q*div(u)*dx #+ gamma*inner(div(u),div(v))*dx
    L = mu*inner(f0, v)*dx + inner(f1,v)*dx
    
    
    
    # Solve
    U = Function(W)
    solve(a == L, U, bcD)
    
    
    # Get sub-functions
    uh, ph = U.split()
    
    
    # error
    
    u = Expression(("pow(x[0],2)*pow(1-x[0],2)*x[1]*(1-x[1])*(1-2*x[1])",
                      "-x[0]*(1-x[0])*(1-2*x[0])*pow(x[1],2)*pow(1-x[1],2)"),degree=10)
    
    p = Expression(("10*(pow(x[0]-1/2,3)*pow(x[1],2)+pow(1-x[0],3)*pow(x[1]-1/2,3))"),degree=6)
    
  

    erroru_H10 = errornorm(u, uh, "H10")
    # erroru_Hdiv0 = errornorm(u,uh,"Hdiv0")
    erroru_Hdiv0 = norm(uh,"Hdiv0")
    return erroru_H10,erroru_Hdiv0
    
N = [2,4,8,16,32,64]
mu = [1,1e-2,1e-4,1e-6]
gamma = [1e-4,1e-2,1,10]

for k in range(0,4,1):
    l = np.size(N)
    s = np.size(mu)
    err = np.zeros([l,s])
    errHdiv = np.zeros([l,s])
    rate = np.zeros([l,s])
    rateHdiv = np.zeros([l,s])
    h = np.zeros(l)
    for i in range(0,l,1):
        h[i] = 1.0/N[i]
        for j in range(0,s,1):
            errs = Stokes(N[i],mu[j],gamma[k])
            err[i,j]= errs[0]
            errHdiv[i,j] = errs[1]
            if (i!=0):
                rate[i,j] = -log(err[i,j]/err[i-1,j])/log(2)
                rateHdiv[i,j] = -log(errHdiv[i,j]/errHdiv[i-1,j])/log(2)
         
    print("\t\t err of H10 of gamma=%.1e"%gamma[k])    
    print("  \tmu \t    1\t\t    1e-2\t     1e-4\t   1e-6 ")
    print("  h\t       err     rate\t err    rate\t   err   rate\t  err  rate")
    i = 0
    print("%.2e"%h[i], " &  %.2e"%err[i,0], " &  ---","  & %.2e"%err[i,1]," & ---"," &   %.2e"%err[i,2]," &---"," &%.2e"%err[i,3],"&---\\\\")
    
    for i in range(1,l,1):
        print("%.2e"%h[i], "  & %.2e"%err[i,0], " & %.2f"%rate[i,0],"  & %.2e"%err[i,1]," & %.2f"%rate[i,1]," & %.2e"%err[i,2],"&%.2f"%rate[i,2],"& %.2e"%err[i,3],"&%.2f\\\\"%rate[i,3])
    
    
    print("\n\n\t\t err of Hdiv0 of gamma=%.1e"%gamma[k])    
    print("  \tmu \t    1\t\t    1e-2\t     1e-4\t   1e-6 ")
    print("  h\t       err     rate\t err    rate\t   err   rate\t  err  rate")
    i = 0
    print("%.2e"%h[i], "  & %.2e"%errHdiv[i,0], "  & ---"," &  %.2e"%errHdiv[i,1],"&  ---","  &  %.2e"%errHdiv[i,2],"& ---","& %.2e"%errHdiv[i,3],"&---\\\\")
    
    for i in range(1,l,1):
        print("%.2e"%h[i], "  & %.2e"%errHdiv[i,0], " & %.2f"%rateHdiv[i,0]," &  %.2e"%errHdiv[i,1],"&  %.2f"%rateHdiv[i,1],"  & %.2e"%errHdiv[i,2],"&%.2f"%rateHdiv[i,2],"& %.2e"%errHdiv[i,3],"&%.2f\\\\"%rateHdiv[i,3])
    
