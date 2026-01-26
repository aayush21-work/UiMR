
import matplotlib
matplotlib.use("Qt5Agg")  # or "QtAgg"
import matplotlib.pyplot as plt 
import numpy as np 

shape=np.loadtxt('shape.dat',dtype=int)
nt,rows,cols=(500,100,100)
nt=nt-1
A=np.loadtxt('solution.dat')
A=A.reshape(nt+1,rows,cols)

plt.figure()
im = plt.imshow(A[0], origin="lower", cmap="inferno")
plt.colorbar()
plt.title("Iteration 0")

for k in range(nt):
    im.set_data(A[k])
    plt.title(f"Iteration {k}")
    plt.pause(0.05)

plt.show()
#print(A)
