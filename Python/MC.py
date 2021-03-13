import numpy as np
import scipy.integrate as integrate
import sys
import matplotlib.pyplot as plt
from scipy import integrate


#Trapezoidal parameters
def Trap(f,n,a,b):
   t = (b-a)/float(n)
   y = 0
   x = a
   for i in range(1,int(n),1):
       x = x+t
       y = y+ f(x)
   y = 0.5*(f(a)+f(b)) +y
   return t*y

#Gauss_Quadrature parameter
def Gauss_quad_lag(n,func,a,b):
    value = 0
    x, w = np.polynomial.legendre.leggauss(int(n))
    for i in range(1,int(n),1):
        value = value+ 0.5*(b-a)*w[i]*func(0.5*(a + b) + 0.5*(b-a)*x[i])
    return value


 #Monte_Carlo integration

def Monte_Carlo(n,func,a,b):
   sub = np.arange(0,n+1,n/10000)
   A = np.zeros(n)
   for i in range(10000):
      ll = int(sub[i])
      ul = int(sub[i+1])
      A[ll:ul] = np.random.uniform(i/10000, (i+1)/10000, ul-ll)
   X = func(a + (b-a)*A)
   MC_integ = ((b-a)/n)*(X.sum())
   return MC_integ


#user defined function
def func(x):
    return 4*x*x*x*np.exp(-x)



partition = 100
a = 0.
b = 20.


if '-limit' in sys.argv:
    p = sys.argv.index('-limit')
    a = int(sys.argv[p+1])
    b = int(sys.argv[p+2])
if '-partition' in sys.argv:
    p = sys.argv.index('-partition')
    partition = int(sys.argv[p+1])

if '-h' in sys.argv or '--help' in sys.argv:
    print ("Usage: %s [-limit] lowerlimit upperlimit [-partition] number " % sys.argv[0])
    print
    sys.exit(1)



    

Trapezoidal = Trap(func,partition,a,b)
gaus_quad = Gauss_quad_lag(partition,func,a,b)
Monte_Carlo_int = Monte_Carlo(partition,func,a,b)
print( f"Trapezoidal Rule a= {a} b = {b} and n = {partition}: ",Trapezoidal)
print (f"Gaussian-Quadrature n ={partition}: ", gaus_quad)
print (f"Monte_Carlo n = {partition}: ", Monte_Carlo_int)
