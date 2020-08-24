from matplotlib import pyplot as plt
from math import *

def ext(x):
    y = -(x+1)*exp(-2*x) +1 
    return y
with open('./UH_hydrogen1.dat','r') as f:
    data = f.read().split('\n')
    data = [i for i in data if i]
    X = [float(i.split()[0]) for i in data]
    Y = [float(i.split()[1]) for i in data]
    plt.scatter(X,Y,c='red')
    plt.plot(X,[ext(i) for i in X],label='ext')
    plt.xlabel('Distance from hydrogen atom')
    plt.ylabel('r$\\times$Vh')
    plt.legend()
    plt.show()


