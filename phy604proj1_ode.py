#PHY607 Project 1 ODE portion

#imports are numpy, matplotlib, sympy

import numpy as np
import matplotlib.pyplot as plt
from sympy import Symbol

#2.0 The differential equation half---

#2.1 Choose ODE method solver
E=0.0301; RK4=0.0302; Both=0.0303
B = eval(input("Please input which integration method - E, RK4, - or Both"))

#2.2 The differential equation of choice is-
#u*d^2/dt^2 r(t) = L^2/(u*r(t)^3)-G*M*u/r(t)^2

#2.21 ICs : Let r(0)=r0, d/dt r(t)|t=0 = v0, the update equations:
r0, v0 = map(float,input("Input initial position, velocity: ").split(","))

#2.22 Precision of diff. eq (steps)
Nd, dt = map(float,input("Input desired step quantity, and step size: ").split(","));Nd=int(Nd)
timee=np.empty(Nd+1); timee[0]=0
r=np.empty(Nd+1);r[0]=r0
rdot=np.empty(Nd+1);rdot[0]=v0

#L=Symbol('L'); u=Symbol('u');
#G=Symbol('G'); M=Symbol('M')
#^having float based issues... don't have more time/not sure how to fix
#input parameters
G, M, m, thetadot = map(float,input("Input Newton's constant, large object mass M, small object m, thetadot: ").split(","))
u=M*m/(M+m); L=u*r0**2*thetadot

if B==E:
  for i in range(Nd):
    timee[i+1]=timee[i]+dt
    rdot[i+1]=float(rdot[i]+dt*(L**2/(u**2*r[i]**3)-G*M/r[i]**2))
    r[i+1]=float(r[i]+dt*rdot[i+1])
    if i % 100 == 0:
      print("r, rdot (aka V) at time t=",timee[i],":",r[i],rdot[i])
    if r[i]<=0:
      print("hit physical barrier-input parameters predict collision in finite time")
if B==RK4:
  def accel(r):
    return (L**2/(u**2*r**3)-G*M/r**2)
  for i in range(Nd):
    #thinking about below was an exercise (for the better or worse)
    a1=dt*rdot[i]
    b1=dt*accel(r[i])
    a2=dt*(rdot[i]+b1/2)
    b2=dt*accel(r[i]+a1/2)
    a3=dt*(rdot[i]+b2/2)
    b3=dt*accel(r[i]+a2/2)
    a4=dt*(rdot[i]+b3)
    b4=dt*accel(r[i]+a3)
    #using a, b to update, oh and split d^2/dt^2 r = d/dt v
    timee[i+1]=timee[i]+dt
    r[i+1]=r[i]+(1/6)*(a1+2*a2+2*a3+a4)
    rdot[i+1]=rdot[i]+(1/6)*(b1+2*b2+2*b3+b4)
    if i % 100 == 0:
      print("r, rdot (aka V) at time t=",timee[i],":",r[i],rdot[i])
    if r[i]<=0:
      print("hit physical barrier-input parameters predict collision in finite time")
#plt.plot(r,timee)
#to get an orbit equation, by Ldot=0, L=u*r**2*thetadot => theta(t)=thetadot*t (for theta0=0)
plt.figure(1)
plt.title(f"radius vs theta mod 2pi, theta0=0, r0={r0},v0={v0}")
plt.plot((thetadot*timee % (2*np.pi)),r)

plt.figure(2)
plt.title(f"Orbit plot position (x,y) versus time, theta0=0, r0={r0},v0={v0}")
plt.plot(r*np.sin(thetadot*timee),r*np.cos(thetadot*timee))
plt.show()

#for i in range(Nd):
#  print(thetadot*timee[i] % (2*np.pi))

#----
#Generate differential equation comparisons



#-------------------------------------

#Analysis
