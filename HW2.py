import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as pt

def condition_number(A):
    _, s, _ = la.svd(A)
    return s[0] / s[-1]

x_normal_list = []
x_qr_list = []
r_normal_list = []
r_qr_list = []
cond_A_list = []

for A, b, x_accurate in zip(A_list, b_list, x_accurate_list):

    # Add code here
    q,r=la.qr(A)
    x_normal = la.solve(A.T@A,A.T@b)
    x_qr = la.solve(r, q.T@b)
    r_normal = la.linalg.norm(x_normal-x_accurate)/la.linalg.norm(x_accurate)
    r_qr = la.linalg.norm(x_qr-x_accurate)/la.linalg.norm(x_accurate)
    cond_A = condition_number(A)

    x_normal_list.append(x_normal)
    x_qr_list.append(x_qr)
    r_normal_list.append(r_normal)
    r_qr_list.append(r_qr)
    cond_A_list.append(cond_A)
#print(A_list)
#print(b_list)
print ("x_accurate:")
print (*x_accurate_list,sep="\n")
print("x_normal:")
print(*x_normal_list,sep="\n")
print("x_qr:")
print(*x_qr_list,sep="\n")
print("r_normal:")
print(*r_normal_list,sep="\n")
print("r_qr:")
print(*r_qr_list,sep="\n")
#It can be seen that QR Factorization generated more accurate solution than the Normal Equations. Normal Equations method is fast but 
# unstable. Cost is n^3+(1/3)n^3. QR factorization is more efficient. Cost is 2n^3. Also QR factorization has smaller relative error. 
# Add plot code here
pt.title("Relative error vs condition number ")
pt.xlabel("$cond_A$")
pt.ylabel("Relative error ")
pt.loglog(cond_A_list, r_normal_list, label=r"Normal Equations")
pt.loglog(cond_A_list, r_qr_list, label="QR factorization")
pt.legend(loc='upper left')
pt.show()
