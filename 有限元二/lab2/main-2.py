#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 10:24:42 2019

@author: chen
"""



from dolfin import *
import numpy as np
from math import *




# Load mesh
def Stokes(N,Elemtype,mu):
    mesh = UnitSquareMesh(N, N)
    
     # Build function space
    if (Elemtype == "MINI"):
        P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
        B = FiniteElement("Bubble", mesh.ufl_cell(), mesh.topology().dim() + 1)
        V = VectorElement(NodalEnrichedElement(P1, B),dim=2)
        Q = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    elif (Elemtype == "P_2-DGP_0"):
        V = VectorElement("Lagrange", mesh.ufl_cell(), 2)
        Q = FiniteElement("DG", mesh.ufl_cell(), 0)    
    elif (Elemtype =="Taylor-Hood3"):
        k=4
        V = VectorElement("Lagrange", mesh.ufl_cell(), k,dim=2)
        Q = FiniteElement("Lagrange", mesh.ufl_cell(), k-1)     
    else:
        V = VectorElement("Lagrange", mesh.ufl_cell(), 4)
        Q = FiniteElement("DG", mesh.ufl_cell(), 3) 
        
    W = FunctionSpace(mesh, V*Q)
    
    
    
    # Boundaries
    def Dirichlt_boundary(x,on_boundary):
        return x[0] < DOLFIN_EPS or x[0] > 1-DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS
    
    # g_D
    g_D = Expression(("0","0"),degree=3)
    
    #Boundaries conditions
    bcD = DirichletBC(W.sub(0), g_D, Dirichlt_boundary)
    
    
    
    #Source term f
    f = Expression(("0","3*pow(x[1],2)-x[1]+1"), degree=2)
    
    
   
    
    # Construct function space
    (u, p) = TrialFunctions(W)
    (v, q) = TestFunctions(W)
    
    # Define stress and strain
    def epsilon(u):
        return 0.5*(grad(u)+grad(u).T)
    
    # Define variational problem
    a = 2*mu*inner(epsilon(u), epsilon(v))*dx - div(v)*p*dx + q*div(u)*dx
    L = inner(f, v)*dx
    
    
    
    # Solve
    U = Function(W)
    solve(a == L, U, bcD)
    
    
    # Get sub-functions
    uh, ph = U.split()


    # error
    
    u = Expression(("0","0"),degree=10)
    p = Expression(("pow(x[1],3)-0.5*pow(x[1],2)+x[1]-7.0/12.0"),degree=10)
    
   
    
    
    erroru_H01 = errornorm(u, uh, "H10")
    return erroru_H01


#main
Elemtype=["MINI", "P_2-DGP_0", "Taylor-Hood3", "P_4-DGP_3"]
#"MINI", "P_2-DGP_0", "Taylor-Hood3", "P_4-DGP_3"
for os in range(0,4,1):
    mu = [1,1e-2,1e-4,1e-6]
    N = [2,4,8,16,32]
    l = np.size(N)
    s = np.size(mu)
    err = np.zeros([l,s])
    rate = np.zeros([l,s])
    h = np.zeros(l)
    for i in range(0,l,1):
        h[i] = 1.0/N[i]
        for j in range(0,s,1):
            err[i,j]=Stokes(N[i],Elemtype[os],mu[j])
            if (i!=0):
                rate[i,j] = -log(err[i,j]/err[i-1,j])/log(2)
                
         
    print("\t\tThe error of",Elemtype[os],"elem")       
    print("  \tmu \t    1\t\t    1e-2\t     1e-4\t   1e-6 ")
    print("  h\t       err     rate\t err    rate\t   err   rate\t  err  rate")
    i = 0
    print("%.2e"%h[i], " &  %.2e"%err[i,0], " &  ---"," &  %.2e"%err[i,1]," & ---","   & %.2e"%err[i,2],"& ---"," & %.2e"%err[i,3],"& ---\\\\")
    
    for i in range(1,l,1):
        print("%.2e"%h[i], " &  %.2e"%err[i,0], " & %.2f"%rate[i,0]," &  %.2e"%err[i,1],"&  %.2f"%rate[i,1],"  & %.2e"%err[i,2],"& %.2f"%rate[i,2],"& %.2e"%err[i,3],"&%.2f\\\\"%rate[i,3])



