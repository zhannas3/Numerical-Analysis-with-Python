import numpy as np
from numpy import linalg as la
G=np.zeros([P.shape[0],P.shape[1],2]) #define mxnx2 matrix 
m=len(P)
g=np.zeros(2)
for i in range(m):
    n=len(P[i])
    for k in range(n):
        if k!=0 and k!=(n-1):#this corresponds to i=0 and i=n where our gradients are different
            g=np.zeros(2)
            for j in obstacles:
                g=g-2*c1*(P[i][k]-j)/(eps+la.norm(P[i][k]-j)**2)**2 #formula for gradient that was derived in Part 1
            g=g+2*c2*(2*P[i][k]-P[i][k+1]-P[i][k-1]) #add the second term of the gradient
            G[i][k]=g #set each value of the matrix to g
print(c1)
print(P)
print(P.shape)
print(obstacles.shape)
print(G)
