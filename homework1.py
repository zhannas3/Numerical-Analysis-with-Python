import numpy as np
import numpy.linalg as la

network_mat=np.zeros((nnodes,nnodes))
for i,d in enumerate(circuit):
    for j,w in d.items():
        network_mat[i][j]=-1.0/w
        network_mat[j][i]=-1.0/w
        network_mat[i][i]+=1.0/w
        network_mat[j][j]+=1.0/w
print(network_mat)

mat=np.array(network_mat)
rhs=np.zeros(nnodes)
for i in fixed_voltages:
    rhs[i]=fixed_voltages[i]
    for j in range(0,nnodes):
        mat[i][j]=0.0
        mat[i][i]=1.0
print(rhs)
print(mat)

voltages=np.linalg.solve(mat,rhs)
print(voltages)
# Feel free to run this. It will show you the circuit we are considering.
# (This changes every time.)
plot_circuit(circuit, fixed_voltages)
# Once you've found the voltages, you may use this version:
plot_circuit(circuit, fixed_voltages, voltages)
