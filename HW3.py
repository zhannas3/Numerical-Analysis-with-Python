import numpy as np
import numpy.linalg as la
def lanczos(L, x0, niter):
    #Implement the body of this function
    n=L.shape[0]
    Q=np.zeros((n,niter))    #niter=k Q has shape nxk
    T=np.zeros((niter,niter))    #T:kxk triangular matrix 
    one=np.ones(n)      #vector of ones
    one=one/la.norm(one) #normalization 
    x0=x0-one*(one @ x0) #remove the vector of ones from initial vector
    a=np.zeros(niter)  #main diagonal elements of T
    b=np.zeros(niter+1)   #subdiagonal elements of T
    q=x0/la.norm(x0) #normalize initial vector
    qn=np.zeros(n)#initialize
    for k in range (niter):
        Q[:,k]=q[:]
        u=L @ Q[:,k]
        a[k]=q[:].T @u
        u=u-b[k]*qn-a[k]*Q[:,k] #updating the value of u
        b[k+1]=la.norm(u)
        qn[:]=Q[:,k]  #setting to Q
        q=u/b[k+1]
        T[k][k]=a[k] #set the diagonal elements 
        T[k-1][k]=b[k] #and subdiagonal 
        T[k][k-1]=b[k]
    return Q, T

def fiedler_ritz(Q,T):
    #Implement the body of this function
    val, vec=la.eig(T)
    ind=np.argsort(val) #sort to find the smallest value
    vec=vec[:,ind]
    fiedlerVec=Q@vec[:,0] #matrix vector product 
    return fiedlerVec

def fiedler(G, k):
    """
    Calculate the fiedler vector of the graph Laplacian matrix
    'G' using 'k' niter of Lanczos algorithm.
    """
    n, m = G.shape

    assert (n == m), "Matrix should be square !!"

    x0 = np.linspace(1, n, num = n)

    ## You should complete this Lanczos function
    Q, T = lanczos(G, x0, k)

    ## You should complete this Fiedler vector computation function
    fiedlerVec = fiedler_ritz(Q,T)

    partitionVec = np.zeros_like(fiedlerVec)
    mfeidler = np.ma.median(fiedlerVec)

    for i in range(n):
        if (fiedlerVec[i] >= mfeidler):
            partitionVec[i] = 1
        else:
            partitionVec[i] = -1

    return partitionVec

points, triangles = readmesh("mesh.1")
plotmesh(points, triangles)
G = mesh2dualgraph(triangles)
partitionVec = fiedler(G, 150)
plotmesh(points, triangles, partitionVec)
