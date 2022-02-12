import numpy as np
import scipy.linalg as la
#define our system of nonlinear equations 
def f(a,b):
    return np.array([a**4+b**4-1,a**3-b])
#I calculated the Jacobian by hands and write it as a function below
def J(a,b):
    return np.array([[4,4],[3,-1]])
#created a function that can give values for n iterations
def broyden(x,y,ff,B,n):
    guesses=np.zeros((5,2))
    B=J(x,y)   #initial value of B is Jacobian
    f=ff(x,y)
    for i in range(n):
        s=la.solve(B,-1*f) #derived from solving for the zero of Taylor series
        x=x+s[0] #calculate next values
        y=y+s[1]
        diff=ff(x,y)-f #(y_k)
        B=B+(np.outer((diff-np.dot(B,s)),s))/(np.dot(s,s)) #B_(k+1)
        f=ff(x,y) #updating
    return x,y

#below I separately find for each iteration
guesses=np.array([broyden(1,1, f, J, 1 ),broyden(1,1, f, J, 2 ),broyden(1,1, f, J, 3 ),broyden(1,1 ,f,J, 4),broyden(1,1, f, J, 5 )])
