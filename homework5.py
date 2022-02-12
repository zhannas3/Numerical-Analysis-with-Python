import numpy as np
n = np.random.randint(5,20)
w1 = 1.0 + np.random.rand()
w2 = 10.0 + np.random.rand()
import matplotlib.pyplot as plt
F_0 = 0.05

def F(t, w):
    return - (w ** 2.0) * F_0 * np.cos(w * t) * np.ones((n,))

k = np.ones(n + 1) * 1000.0
m = np.ones(n) * 10.0
k[-1] = 0.0
K = np.diag(-k[0:-1] - k[1:]) + np.diag(k[1:-1], 1) + np.diag(k[1:-1], -1)
def f(y, t, w):
    m = int(len(y) / 2)
    x = y[:m]
    v = y[m:]
    z = (K @ x) / 10.0 + F(t, w)
    ynew = np.hstack((v, z)).T
    return ynew
x1=np.zeros(n)
x2=np.zeros(n)

def rk4(t, h, n, IV,y,w):
    y[0]=IV
    y = np.zeros(2*n,float);
    t=np.linspace(0,4,400)
    for i in range(400):
        k1 = f(y,t[i],w)
        k2 = f(y+0.5*k1*h, t[i]+0.5*h, w)
        k3 = f(y+0.5*k2*h, t[i]+0.5*h, w)
        k4 = f(y+k3*h, t[i]+h, w)
        
        y = y+(k1+2*k2+2*k3+k4)*h/6.0
        
    return y[:n]      #return only first n values from Y: displacement 

y_val = np.zeros(2*n)
t=np.linspace(0,4,400)
h = 4/399
print(K.shape)
x1=rk4(t,h,n,0,y_val,w1)
x1=np.array(x1)
print(x1)
print(y_val.shape)
x2=rk4(t,h,n,0,y_val,w2)
x2=np.array(x2)
print(x2)
def plot(x):
    plt.figure()
    xplot = np.insert(x,0,0)
    plt.plot(xplot,np.arange(0,n+1),'o-')
    plt.xlim([-1,1])
    plt.title('$\omega$')
    plt.xlabel('Displacement')
    plt.ylabel('Floors')
    plt.show()
print(n,w1,w2)
