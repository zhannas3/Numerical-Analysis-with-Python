import numpy as np
import matplotlib.pyplot as plt

def not_allowed(*args, **kwargs):
    raise RuntimeError("You called an illegal function.")

import scipy.optimize as opt
for attr in dir(opt):
    setattr(opt, attr, not_allowed)

import scipy.linalg as sla
for attr in dir(sla):
    setattr(sla, attr, not_allowed)

def secular(v,d,beta):
    def f(x):
        return 1 + beta * np.sum((v**2)/(d-x))
    return f

def plot_func(f, x_min,x_max, y_min = -6, y_max = 6, approx_roots = []):
    plt.figure()
    xs = np.linspace(x_min,x_max,1000)
    ys = np.array([f(x) for x in xs])

    # https://stackoverflow.com/q/10377593
    pos = np.where(np.abs(np.diff(ys)) >= 10)[0]+1
    xs = np.insert(xs, pos, np.nan)
    ys = np.insert(ys, pos, np.nan)

    plt.plot([xs[0],xs[-1]], [0,0])
    plt.plot(xs,ys)
    plt.plot(approx_roots, [f(a) for a in approx_roots], "o")
    plt.grid()
    plt.ylim(y_min,y_max)
    plt.show()

def interp(c1,c2,c3,dl,dr,x):
    interpol = c3 + c1/(dl - x) + c2/(dr - x)
    return interpol

def plot_interp(f, c1,c2,c3,dl,dr, x_min,x_max, y_min = -6, y_max = 6):
    plt.figure()
    xs = np.linspace(x_min,x_max,1000)
    ys = np.array([f(x) for x in xs])
    interp_ys = interp(c1,c2,c3,dl,dr,xs)

    # https://stackoverflow.com/q/10377593
    pos = np.where(np.abs(np.diff(ys)) >= 10)[0]+1
    xs = np.insert(xs, pos, np.nan)
    ys = np.insert(ys, pos, np.nan)
    interp_ys = np.insert(interp_ys, pos, np.nan)

    plt.plot([xs[0],xs[-1]], [0,0], label="x-axis")
    plt.plot(xs,ys,label="original function")
    plt.plot(xs,interp_ys,'--', label="approximation")
    plt.ylim(y_min,y_max)
    plt.legend()
    plt.grid()
    plt.show()
#def secular_solve(v,d,beta,idx):
    #return [(d[idx] + d[idx+1])/2, ]
def secular_solve(v,d,beta,idx):
    lambdal=0.0 #our previous estimate
    guesses=[]
    lambda_init=(d[idx]+d[idx+1])/2 #initial guess
    guesses.append(lambda_init)
    n=len(v)
    while not np.isclose(lambda_init,lambdal,rtol=1e-10, atol=1e-14) :
        psi1=0;psi2=0;psi1d=0;psi2d=0
        for j in range(0,idx+1): #j=1 to i
            diff = d[j] - lambda_init
            psi1+=np.power(v[j],2)/np.power(diff,2)
            psi1d+=np.power(v[j],2)/diff
        for j in range(idx+1,n): #j=i+1 to n
            diff = d[j] - lambda_init
            psi2+=np.power(v[j],2)/np.power(diff,2)
            psi2d+=np.power(v[j],2)/diff
        c1=beta*psi1*(d[idx]-lambda_init)**2 #derived formula for c1
        c2=beta*psi2*(d[idx+1]-lambda_init)**2 #derived formula for c2
        c1hat=beta*psi1d-beta*psi1*(d[idx]-lambda_init) #we compute c1hat and c2hat as the following
        c2hat=beta*psi2d-beta*psi2*(d[idx+1]-lambda_init)
        c3hat=1
        c3=1+c1hat+c2hat         #below we compute two zeros of g

        a=c3                           #simplify the coefficients of quadratic equation
        c=c1*d[idx+1]+c2*d[idx]+c3*d[idx]*d[idx+1]
        b=-d[idx]*c3-d[idx+1]*c3-c2-c1
        root1=(-b+(b**2-4*a*c)**0.5)/(2*a) #we find two roots of quadratic equation
        root2=(-b-(b**2-4*a*c)**0.5)/(2*a)
        lambdal=lambda_init
        if root2>d[idx] and root2<d[idx+1]:#check the location of our root
            lambda_init=root2
        else:
            lambda_init=root1
        guesses.append(lambda_init)
    return guesses

beta = 1
v = np.array([0.6,]*4)
d = np.arange(1,5)

f = secular(v,d,beta)
rts = secular_solve(v,d,beta,1)

#For the first iteration of the above call:
c1 = 0.4
c2 = 0.4
c3 = 1.0

#Plot function and correct interpolant for 1st Newton's iteration (using above c's)
plot_interp(f,0.4,0.4,1.0,d[1],d[2], 0,5,y_min = -3, y_max = 3)

###Below are parameters to test your secular_solve function to check for correctness
###uncomment to test

beta_e = 1
v_e = np.array([0.6,]*4)
d_e = np.arange(1,5)
i_e = 1
approx_e = secular_solve(v_e, d_e, beta_e, i_e)

beta_h = 1
v_h = np.array([0.001,]*4)
d_h = np.arange(1,5)
i_h = 2
approx_h = secular_solve(v_h, d_h, beta_h, i_h)

beta_c = 0.8
approx_beta = secular_solve(v_h,d_h,beta_c,i_h)

beta_o = 0.9
v_o = np.array([0.2,0.3,0.004,0.1])
d_o = np.array([-0.2,0.3,1,4])
i_o = 1
approx_o = secular_solve(v_o, d_o, beta_o, i_o)
