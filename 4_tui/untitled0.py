#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
from sympy import Symbol, solve, lambdify, Matrix, Eq
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def f(x, t, k1, k_1, k2, k_2, k3):
    print x
    print 't'
    print t
    print k1
    print k_1
    print k2
    print k_2
    print k3
    xx=x[0]
    y=x[1]
    f1=k1*(1-xx-y)-k_1*xx-k3*xx*y
    f2=k2*(1-xx-y)**2-k_2*y**2-k3*xx*y
    return [f1, f2]
  

k1 = Symbol("k1")
k2 = Symbol("k2")
k3 = Symbol("k3")
k_1 = Symbol("k_1")
k_2 = Symbol("k_2")
x = Symbol("x")
y = Symbol("y")

x_form_y=k1*(1-y)/(k1+k_1+k3*y)
k2_form_y=(k_2*(y**2)*(k1+k_1+k3*y)**2+(k1+k_1+k3*y)*k1*k3*y*(1-y))/((1-y)**2*(k_1+k3*y)**2)
x_fr_y_fun=lambdify((y, k1, k_1, k3), x_form_y)
k2_fr_y_fun=lambdify((y, k1, k_1, k_2, k3), k2_form_y)

y_all = np.linspace(0,0.999,3000)
plt.figure(1)
plt.grid(True)
plt.xlabel('k2')
plt.ylabel('x, y')
plt.xlim((0,8))
plt.ylim((0,1))
x_all=[]
k2_all=[]
for one_y in y_all:
    x_all.append(x_fr_y_fun(one_y, 1, 0.04, 10))
    k2_all.append(k2_fr_y_fun(one_y, 1, 0.04, 0.02, 10))
  
plt.plot(k2_all, x_all, color='r') 
plt.plot(k2_all, y_all, color='b') 
plt.show()


eq1 = k1*(1 - x - y) - k_1*x - k3*x*y
eq2 = k2*(1 - x - y)**2 - k_2*y**2 - k3*x*y
A = Matrix([eq1, eq2])
jacA = A.jacobian(Matrix([x, y]))  
detA = jacA.det()
det_fun = lambdify( (x, y, k1, k_1, k2, k_2, k3), detA)

detA_all=[]
x_bif=[]
y_bif=[]
k2_bif=[]
for i in range(1000):
    detA_all.append(det_fun(x_all[i], y_all[i], 1, 0.004, k2_all[i], 0.02, 10))
    if i>0 and detA_all[i-1]*detA_all[i]<0:
        k2_bif.append(k2_all[i])
        x_bif.append(x_all[i])
        y_bif.append(y_all[i])
    
print k2_bif

y_in_bif = []
for i in y_all:
    k2_res=k2_fr_y_fun(i , 1, 0.04, 0.02, 10)
    if abs(k2_res-k2_bif[0])<0.0001 or abs(k2_res-k2_bif[1])<0.0001:
        y_in_bif.append(i)
        
k2_all1 = []
k2_all2 = []
z = Symbol("z")
z_x_y=1-x-y
k2_form_y_k1=(k1+z)*(2*k_2*y*(k1+z)**2 +(k1+k_1)*k1*k3+k1*k3**2*y**2)/(2*z*(1-y)*(k1*k3*(1-y)-z*(k1+z)))
k2_form_y_k1=k2_form_y_k1.subs(z, z_x_y)
k2_fr_y_k1_fun=lambdify((y, k1, k_1, k3), k2_form_y_k1)

k1_all= np.linspace(0,2,3000)
for i in k1_all:
    k2_all1.append(k2_fr_y_fun(y_in_bif[0], i, 0.04, 0.02, 10))
    k2_all2.append(k2_fr_y_fun(y_in_bif[1], i, 0.04, 0.02, 10))
 
plt.figure(2)
plt.grid(True)
plt.xlabel('k2')
plt.ylabel('k1')
plt.xlim((0,5))
plt.ylim((0,2))
plt.plot(k2_all1, k1_all, color='b', label='1')
plt.plot(k2_all2, k1_all, color='r', label='2')
plt.show()

t = np.linspace(0.01, 1, 1000)
y0=[x_bif[0], y_bif[0]]
#y0=[0, 0]
#sol = odeint(f, y0, t, args=(1, 0.04, k2_bif[0], 0.02, 10))
sol = odeint(f, y0, t, args=(1, 0.04, 2, 0.02, 10))

print k2_bif[0]

plt.figure(3)
plt.grid(True)
plt.xlabel('t')
plt.ylabel('x, y')
plt.plot(t, sol[:, 0], color='b')
plt.plot(t, sol[:, 1], color='g')
plt.show()

plt.figure(4)
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')
plt.plot(sol[:, 0], sol[:, 1], color='r')
plt.show()

#==============================================================================
# 
# z = Symbol("z")
# z_x_y=1-x-y
# k2_form_y_k1=(k1+z)*(2*k_2*y*(k1+z)**2 +(k1+k_1)*k1*k3+k1*k3**2*y**2)/(2*z*(1-y)*(k1*k3*(1-y)-z*(k1+z)))
# k2_form_y_k1=k2_form_y_k1.subs(z, z_x_y)
# k2_fr_y_k1_fun=lambdify((y, k1, k_1, k3), k2_form_y_k1)
# plt.figure(2)
# plt.xlabel('k2')
# plt.ylabel('k1')
# plt.xlim((0,5))
# plt.ylim((0,2))
# k11_all=[]
# k12_all=[]
# for i in range(1000):
#     my_eq=Eq(k2_form_y_k1, k2_all[i])
#     tmp= solve(my_eq, k1)
#     print tmp
#     #k11_all.append(tmp[0].evalf(subs={y: y_all[i], k_1: 0.04, k_2:0.02, k3:10}))
#     #k12_all.append(tmp[1].evalf(subs={y: y_all[i], k_1: 0.04, k_2:0.02, k3:10}))
# ##plt.plot(k2_all, x_all, color='r') 
# #plt.plot(k2_all, y_all, color='b') 
# #plt.show()
#     
#==============================================================================
