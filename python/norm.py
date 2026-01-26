import numpy as np 

file=np.loadtxt('bounds1.dat')
print(file)
print("--------------------")
nx,ny=file.shape

for i in range(0,nx-1):
    tempx = (file[i][1]+file[i+1][0])/2
    file[i][1]=tempx
    file[i+1][0]=tempx

for j in range(0,ny-1):
    tempy = (file[j][2]+file[j][3])/2
    file[j][2]=tempy
    file[j][3]=tempy
print(file)
