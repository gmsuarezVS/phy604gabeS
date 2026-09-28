#PHY607 Project 1

#imports are numpy, matplotlib, sympy

import numpy as np
import matplotlib.pyplot as plt
from sympy import Symbol
from scipy.special import beta as scipy_beta
from matplotlib.ticker import LinearLocator
#1.0 The integral half ---

#1.1 choose integral type
R=0.0501; T=0.0502; S=0.0503;
A = eval(input("Please input which integration method - R, T, S: "))
#Note: no input for function to integrate-I've chosen it, but could be done

#1.2 Here we keep s symoblic - so we can plot in s, t after we integrate out z
def fi(x,s,t):
  if x <= 0: x = 0.0001
  if x >= 1: x = 1- 0.0001
  return x**(-s-1)*(1-x)**(-t-1)

#1.3 For concrete case by case, define some s1, t1
s1, t1 = map(float,input("Please input your Mandelstam variables (expect convergence, generically, only for s,t < 0): ").split(","))
print(s1,t1)
def gi(z,s1,t1):
  if z <= 0: x = 0.0001
  if z >= 1: x = 1 - 0.0001
  return z**(-s1-1)*(1-z)**(-t1-1)

#1.4 Fill out values of fi, gi

#1.41 input desired precision on integrand range
N = int(input("Input how many points you want to sample between [0,1]: "))
xi=np.empty(N); xi[0]=0
for i in range(N):
  xi[i]=0+i/N


#1.42 for later - sample over s, t space
Ns, Nt = map(int,input("Input how many points you want to sample s,t over: ").split(","))
Rs1, Rs2, Rt1, Rt2 = map(float,input("Input range of sampling s,t (if you want -1,1 for both, input '-1, 1, -1, 1'): ").split(","))
si=np.empty(Ns+1)
for i in range(Ns+1):
  si[i]=Rs1+(Rs2-Rs1)*i/Ns

ti=np.empty(Nt+1) ;
for i in range(Nt+1):
  ti[i]=Rt1+(Rt2-Rt1)*i/Nt

#1.5 Time to do symbolic integration over fi, gi
s =Symbol('s')
t =Symbol('t')

#1.51 define Riemann integration
if A==R:
    Rf=0; Rg=0
    #print("si",si,"ti",ti)
    for i in range(N):
      Rf = Rf+fi(i/N,s,t)*(1/N)
      Rg = Rg+gi(i/N,s1,t1)*(1/N)
    #print("Rf:",Rf,"Rg:",Rg)
#1.52 define Trapezoidal integration
elif A==T:
  Rf=0; Rg=0
  for i in range(N+1):
    Rf = Rf+(fi((i+1)/N,s,t)+fi(i/N,s,t))*(1/2)*(1/N)
    Rg = Rg+(gi((i+1)/N,s1,t1)+gi(i/N,s1,t1))*(1/2)*(1/N)
  #print("Tf:",Rf,"Tg:",Rg)
#1.53 define Simpson's integration
elif A==S:
  Rf=fi(0,s,t)*1/(N*3)+fi(1,s,t)*1/(N*3); Rg=gi(0,s1,t1)*1/(N*3)+gi(1,s1,t1)*1/(N*3)
  for i in range(1,N):
    if i % 2 == 0:
      Rf=Rf+2*fi(i/N,s,t)*1/(N*3)
    elif i % 1==0:
      Rf=Rf+4*fi(i/N,s,t)*1/(N*3)
  for i in range(1,N):
    if i % 2 == 0:
      Rg=Rg+2*gi(i/N,si,ti)*1/(N*3)
    elif i % 1==0:
      Rg=Rg+4*gi[i/N,si,ti]*1/(N*3)
  #print("Sf:",Rf,"Sg:",Rg)
#1.55 to show typo
else:
  print("You put A=",A,",which is not a valid choice... ")

#----
#1.6 Generate integration comparison printouts
print("This compares the computed value of the integral, for some, s,t. EB:",scipy_beta(-s1,-t1),"computed:",Rg)
RRf=np.empty((Ns+1,Nt+1)); RRB=np.empty((Ns+1,Nt+1))
for i in range(Ns+1):
  Rfs=Rf.subs(s,si[i])
  for j in range(Nt+1):
    RRf[i,j]=float(Rfs.subs(t,ti[j]))
    RRB[i,j]=float(scipy_beta(-si[i],-ti[j]))
#print("Computed:",RRf,"Scipy:",RRB,"Differences:",RRf-RRB)

#1.7 Generating 3D Plot comparisons over s,t
#1.71 comparison grids per https://numpy.org/doc/2.5/reference/generated/numpy.meshgrid.html
sgrid, tgrid = np.meshgrid(si,ti,indexing='ij')

#plot 1 per https://matplotlib.org/stable/gallery/mplot3d/surface3d.html
fig1, ax1 = plt.subplots(subplot_kw={"projection": "3d"})
# Make data.
ax1.set_title(f"Numerical integration method B(s,t)")
ax1.set_xlabel('s') ;ax1.set_ylabel('t');
# Plot the surface.
surf = ax1.plot_surface(sgrid, tgrid, RRf, cmap="coolwarm",
                       linewidth=0, antialiased=False)
# Add a color bar which maps values to colors.
fig1.colorbar(surf, shrink=0.5, aspect=5)


#----

#plot 2 
fig2, ax2 = plt.subplots(subplot_kw={"projection": "3d"})
# Make data.
ax2.set_title(f"SciPy integration method B(s,t)")
ax2.set_xlabel('s') ;ax2.set_ylabel('t');
# Plot the surface.
surf = ax2.plot_surface(sgrid, tgrid, RRB, cmap="coolwarm",
                       linewidth=0, antialiased=False)
# Add a color bar which maps values to colors.
fig2.colorbar(surf, shrink=0.5, aspect=5)


#----

#plot 3 
fig3, ax3 = plt.subplots(subplot_kw={"projection": "3d"})
# Make data.
ax3.set_title(f"Numerical % Difference from SciPy B(s,t)")
ax3.set_xlabel('s') ;ax3.set_ylabel('t');
# Plot the surface.
surf = ax3.plot_surface(sgrid, tgrid, 100*(RRB-RRf)/RRB, cmap="coolwarm",
                       linewidth=0, antialiased=False)
# Add a color bar which maps values to colors.
fig3.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

# analysis for integration
