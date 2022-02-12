import numpy as np
import matplotlib.pyplot as pt

def f(x):
    return np.sin(k * x)
def ft(x):
    return -(k**2)*np.sin(k*x)

h_values = []
abs_err_values = []
truncation_err_values=[]
rounding_err_values=[]
for n in n_values:
    x = np.linspace(0, 1, n)
    eps = np.finfo(np.float64).eps
    fx = f(x)
    ftx = ft(x)
    h=1/(n-1)
    h_values.append(h)
    dfx_app = (f(x+h)-2*f(x)+f(x-h)) / h**2
    abs_err_values.append(np.max(np.abs(dfx_app-ftx)))
    truncation_err_values.append(((k**4)*(h)**2)/12)
    rounding_err_values.append(4*eps*k/(h)**2)
abs_err_values = np.array(abs_err_values,dtype=np.float64)
truncation_err_values=np.array(truncation_err_values)
rounding_err_values=np.array(rounding_err_values)

print(abs_err_values)
print(truncation_err_values)
print(rounding_err_values)

# ------------ plotting code below, no need to change
pt.xlabel(r"$h$")
pt.loglog(h_values, truncation_err_values + rounding_err_values, label=r"Predicted Error Bound")
pt.loglog(h_values, abs_err_values, label="Computed Absolute Error")
pt.legend(loc='lower right')
pt.show()
