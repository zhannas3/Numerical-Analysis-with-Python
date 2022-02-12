import numpy as np
n = np.random.randint(100,150)
r = np.random.randint(5,15)

#construct rank-deficient A
A = np.random.rand(n,n)
u,s,vt = np.linalg.svd(A)
s1 = np.array([10**(-i) for i in range(r)])
s[-r:] = s1
A = u @ np.diag(s) @ vt

def apply_householder(Q, R, col):
    #compute Householder vector
    v = np.zeros(R[col:, col].shape)
    v[0] += np.sign(R[col, col])*np.linalg.norm(R[col:, col])
    v += R[col:,col]

    #apply Householder vector
    k = 2./np.inner(v,v)
    R[col:,col:] -= k * np.outer(v, v @ R[col:,col:])
    Q[:,col:] -= k * np.outer(Q[:,col:] @ v, v)

def not_allowed(*args, **kwargs):
    raise RuntimeError("You called an illegal function.")

import scipy.linalg as sla
for attr in dir(sla):
    setattr(sla, attr, not_allowed)
attrs = dir(np.linalg)
attrs.remove('norm')
for attr in attrs:
    setattr(np.linalg, attr, not_allowed)
R=A.copy()
Q=np.eye(A.shape[0])
P=np.arange(A.shape[0],dtype=np.int64)
for i in range(A.shape[0]-1):
    cnorm=np.linalg.norm(R[i:,i:], ord=2, axis=0)
    imax=np.argmax(cnorm)
    P[[i+imax,i]]=P[[i,i+imax]]
    R[:,[i+imax,i]]=R[:,[i,i+imax]]
    H=apply_householder(Q,R,i)
