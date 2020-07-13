#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 21:32:20 2019

@author: chen
"""



from dolfin import *
import numpy as np
from math import *
import os, sys
import matplotlib.pyplot as plt
import sympy as sp






def Stokes(N,Elemtype):
    #exact solution
    
    
    x,y = sp.symbols("x[0],x[1]")
    u1 = sp.exp(x)*sp.cos(sp.pi*y)
    u2 = -1/sp.pi*sp.exp(x)*sp.sin(sp.pi*y)
    ps = (x-0.5)*(y-0.5)
    f1 = -sp.diff(u1,x,2) - sp.diff(u1,y,2) + sp.diff(ps,x,1)
    f2 = -sp.diff(u2,x,2) - sp.diff(u2,y,2) + sp.diff(ps,y,1)
    g1 = 2*sp.diff(u1,x,1)-ps
    g2 = sp.diff(u1,y,1)+sp.diff(u2,x,1)
    # transform the expressions to ccode style
    u1 = sp.printing.ccode(u1)
    u2 = sp.printing.ccode(u2)
    ps = sp.printing.ccode(ps)
    f1 = sp.printing.ccode(f1)
    f2 = sp.printing.ccode(f2)
    g1 = sp.printing.ccode(g1)
    g2 = sp.printing.ccode(g2)
    
    exactu = Expression((u1,u2),degree=10)
    exactp = Expression(ps,degree=10)
    
 
    
    # Load mesh
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
        k = 2
        V = VectorElement("Lagrange", mesh.ufl_cell(), k)
        Q = FiniteElement("Lagrange", mesh.ufl_cell(), k-1)     
    elif (Elemtype =="P_4-DGP_3"):
        V = VectorElement("Lagrange", mesh.ufl_cell(), 5)
        Q = FiniteElement("DG", mesh.ufl_cell(), 4) 
    elif (Elemtype =="P_1-P_1"):
        V = VectorElement("Lagrange", mesh.ufl_cell(), 1)
        Q = FiniteElement("Lagrange", mesh.ufl_cell(), 1) 
         
    
    W = FunctionSpace(mesh, V*Q)
    
    
    # Boundaries
    def Dirichlt_boundary(x,on_boundary):
        return x[0] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS
    
    g_D =Expression((u1,u2),degree=10)
   
    #Boundaries conditions
    bcD = DirichletBC(W.sub(0), g_D, Dirichlt_boundary)
    
    #Neumann g_N

    g_N = Expression((g1,g2), degree=5)
   
    #Source term f
   
    f = Expression((f1,f2), degree=10)
    
    
    # Construct function space
    (u, p) = TrialFunctions(W)
    (v, q) = TestFunctions(W)
    
    # Define stress and strain
    def epsilon(u):
        return 0.5*(grad(u)+grad(u).T)
    
    # Define variational problem
    a = 2*inner(epsilon(u), epsilon(v))*dx - div(v)*p*dx + q*div(u)*dx
    L = inner(f, v)*dx + inner(g_N,v)*ds
    
    
    
    # Solve
    U = Function(W)
    solve(a == L, U, bcD)
    
    
    # Get sub-functions
    u, p = U.split()
    
  
    
    # error
    
    
    
    erroru_L2 = errornorm(exactu, u, "L2")
    erroru_H1 = errornorm(exactu, u, "H1")
    errorp_L2 = errornorm(exactp, p, "L2")
    
    return erroru_L2,erroru_H1,errorp_L2


#main
#Choose elem type
Elemtype=["MINI", "P_2-DGP_0", "Taylor-Hood3", "P_4-DGP_3","P_1-P_1"]
#"MINI", "P_2-DGP_0", "Taylor-Hood3", "P_4-DGP_3","P_1-P_1"
for s in range (0,5,1):
    N = [2,4,8,16,32,64]
    l = np.size(N)
    
    Errorp = np.zeros([l])
    Erroru = np.zeros([l])
    Erroru1 = np.zeros([l])
    
    Ratep = np.zeros([l])
    Rateu = np.zeros([l])
    Rateu1 = np.zeros([l])
    h = np.zeros([l])
    
    for i in range(0,l,1):
        Error=Stokes(N[i],Elemtype[s])
        Erroru[i] = Error[0]
        Erroru1[i] = Error[1]
        Errorp[i] = Error[2]
        if (i!=0):
            Ratep[i] = -log(Errorp[i]/Errorp[i-1])/log(2)
            Rateu[i] = -log(Erroru[i]/Erroru[i-1])/log(2)
            Rateu1[i] = -log(Erroru1[i]/Erroru1[i-1])/log(2)
        h[i] = 1/N[i]
        
    print("\t\tThe error of",Elemtype[s],"elem")       
    print("h          error-u-L2 rate   error-u-H1 rate   error-p-L2 rate")
    i = 0
    print("%.2e"%h[i], " &  %.2e"%Erroru[i], " & ---","  & %.2e"%Erroru1[i]," & ---"," &  %.2e"%Errorp[i],"& ---\\\\")
    for i in range(1,l,1):
        print("%.2e"%h[i], " &  %.2e"%Erroru[i], "&% .2f"%Rateu[i]," &  %.2e"%Erroru1[i],"&% .2f"%Rateu1[i]," &  %.2e"%Errorp[i]," &%.2f\\\\"%Ratep[i])
    
    




