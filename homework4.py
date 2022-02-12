import scipy.optimize as sopt
import numpy as np
import matplotlib.pyplot as pt
x = init_path
final_path=np.zeros(np.shape(x))
iters=[10,25,50,100,200]
#Draw the canvas
pt.figure(figsize=(8,8))
pt.plot(x[:,0], x[:,1], label="Initial")
for i in range(400):
    c1=1000/(100+i)
    c2=1
    s=-dobj(x,c1,c2)
    def f1d(p):
        return (obj(x-p*dobj(x,c1,c2),c1,c2))     
    alpha=sopt.golden(f1d)        
    x=x+alpha*s  
    final_path=x
    if i in iters:
        pt.plot(x[:,0],x[:,1],label=f" the path after {i}th iterations")
pt.plot(final_path[:,0],final_path[:,1], label="final path")

pt.plot(obstacles[:,0], obstacles[:,1], "o", markersize=5, label="Obstacles")
pt.legend(loc="best")
#print(steep(n=400))
#print(x)
#print(obstacles)
