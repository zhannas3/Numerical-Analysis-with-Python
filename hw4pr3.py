import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

def plot_func():
    (a, b) = interval
    xs = np.linspace(a, b, 10_000)
    plt.plot(xs, f(xs))

def plot_interpolant(chunks):
    nodes = []
    interpolant = []
    nodes_biunit = np.linspace(-1, 1, 20)

    basis = [np.ones(20), nodes_biunit]
    for i in range(2, nnodes):
        basis.append(2*nodes_biunit*basis[-1] - basis[-2])
    vdm = np.array(basis).T

    nodes_01 = np.linspace(0, 1, 20)
    for a, b, coeffs in chunks:
        nodes.extend(a+(b-a)*nodes_01)
        interpolant.extend(vdm @ coeffs)

    plt.plot(nodes, interpolant, "o-")
    plt.plot(nodes, f(np.array(nodes)), "o-")
a, b=interval

def adaptive_interp(interval):
    # Add code here...
    a,b=interval
    k=[]
    for i in range(nnodes):#nnodes=n
        k.append(float(i))
    i=np.array([k])
    xk=np.cos((2*(i+1)-1)*np.pi/(2*nnodes)) #formula of Chebyshev nodes for interval (-1,1)
    x=xk*(b-a)/2+(b+a)/2 #for arbitrary [a,b] interval
    T=np.cos(i*np.arccos(xk).reshape(-1,1)) #define Chebyshev polynomial
    l=[]
    for i in np.nditer(x):
        l.append(f(i))
    xleft=np.array(l)
    alpha=la.solve(T, xleft)
    error=np.sqrt(np.power(alpha[-2],2)+np.power(alpha[-1],2))/np.sqrt(np.sum(alpha**2)) #error formula
    if error<tol:
        return [(a,b,alpha)]
    else: #else: we split interval
        return adaptive_interp((a,(b+a)/2))+adaptive_interp(((b+a)/2, b))
    # Making this function recursively call itself may be helpful.
    #pass

interpolant_chunks = adaptive_interp(interval)
print(f"{len(interpolant_chunks)} chunks produced")


# You may use this to see the functions you are working on.
plot_func()

# Once you have interpolant_chunks computed, you may enable this to
# see your interpolant visualized compared to the original function.
plot_interpolant(interpolant_chunks)
