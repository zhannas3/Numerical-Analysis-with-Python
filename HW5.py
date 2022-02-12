import numpy as np
import matplotlib.pyplot as plt
import math as mt
def f(x):
    return 4.0/(1.0+x**2)
#midpoint
def mid(a,b,n): 
    h = (b-a) / float(n)
    mid = 0
    for i in range(n):        
        mid += f((a + h/2.0) + i*h)
    mid *= h
    return mid
intgmid =[]
emid=[]
for i in ks:   
    m = mid(0, 1, i)
    ea=abs((mid(0,1,i)-np.pi) /np.pi)    
    emid.append(ea) 
emid=np.array(emid)
print(emid)
#trapezoid 
def trap(a,b,n):
    h = (b-a) / float(n)
    intgr = 0.5 * h * (f(a) + f(b))
    for i in range(1, int(n)):
        intgr = intgr + h * f(a+i * h)
    return intgr
etrap=[]
for i in ks:    
    t = trap(0, 1, i)
    ea1=abs((trap(0,1,i)-np.pi) /np.pi)    
    etrap.append(ea1)         
etrap=np.array(etrap)
print(etrap)
#simpson can be defined as the linear combination of midpoint and trapezoid
def simp(a,b,n):
    o=mid(a,b,n)
    p=trap(a,b,n)
    return 2/3*o+1/3*p
esimp=[]
for i in ks:    
    ss = simp(0, 1, i)
    ea3=abs((simp(0,1,i)-np.pi) /np.pi)    
    esimp.append(ea3)      
esimp=np.array(esimp)
print(esimp)
print(nevals_gauss)
print(ks)
print(len(nevals_gauss))
#gauss-legendre
def gaussxc(n):
    # Initial approximation to roots of the Legendre polynomial
    a = np.linspace(3,4*n-1,n)/(4*n+2)
    x = np.cos(mt.pi*a+1/(8*n*n*np.tan(a)))
    # Find roots of Legendre polynomial using Newton's method
    e = 1e-16
    delta = 1.0
    while delta>e:
        p0 = np.ones(n,float)
        p1 = np.copy(x)
        for k in range(1,n):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (n+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x =x- dx
        delta = max(abs(dx))
    # Calculate the coefficients
    c = 2*(n+1)*(n+1)/(n*n*(1-x*x)*dp*dp)
    return x,c
#returns the value of the integral with n number of points(n-point Gauss-Legendre Formula)
def gauss(n,a,b):
    x,c = gaussxc(n)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*c
#define a function for error
def gausserror(n):
    a=0
    b=1
    xd,dxd = gaussxc(n)
    x = 0.5*(b-a)*xd + 0.5*(b+a)
    dx = 0.5*(b-a)*dxd   
    egauss=[]
    s =sum(dx*f(x))
    ea3=abs((s-np.pi) /np.pi)             
    return ea3
egauss=np.array([gausserror(nevals_gauss[0]),  gausserror(nevals_gauss[1]),  gausserror(nevals_gauss[2]),  gausserror(nevals_gauss[3]),
                 gausserror(nevals_gauss[4]), gausserror(nevals_gauss[5]), gausserror(nevals_gauss[6]),
                gausserror(nevals_gauss[7]), gausserror(nevals_gauss[8]), gausserror(nevals_gauss[9]), gausserror(nevals_gauss[10])])
print(egauss)
