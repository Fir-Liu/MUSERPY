# -*- coding: utf-8 -*-
"""
Created on Fri Apr 07 09:17:28 2017

@author: super
"""
s1 = 72
s2 = 85

r = 'Xiaoming has improved his gpa by percent %2.1f%%'
p = (s2-s1)/s1*100
print(r % p)

#%%
L = [
    ['Apple', 'Google', 'Microsoft'],
    ['Java', 'Python', 'Ruby', 'PHP'],
    ['Adam', 'Bart', 'Lisa']
]
t = ('a', 'b', ['A', 'B'])
#%%
#height = input('height(m):')
#weight = input('weight(kg):')
#bmi = float(weight)/(float(height)**2)
#
#if bmi < 18.5:
#    sit = 'skinny'
#elif bmi <= 25:
#    sit = 'normal'
#elif bmi <= 28:
#    sit = 'less fat'
#elif bmi <= 32:
#    sit = 'fat'
#else:
#    sit = 'overfat'
#   
#r = 'Your BMI is %2.1f, and you are %s'
#print(r % (bmi, sit))

import matplotlib.pyplot as plt
import numpy as np

height = input('height(m):')
myheight = float(height)
weightgrp = np.linspace(50,69,20,endpoint=True)
bmigrp = weightgrp/myheight
#plt.plot(weightgrp,bmigrp,'o')
fig, ax = plt.subplots()
ax.plot(weightgrp,bmigrp,'o')
ax.grid(True)

#%%
L = ['Bart', 'Lisa', 'Adam']
for i in range(len(L)):
    print('hello %s' % L[i])
#%%
n1 = 255
n2 = 1000
print(hex(n2))
#%%
#import math

def quadratic(a, b, c):
    if a!=0:
        r1 = (-b + (b**2-4*a*c)**0.5)/(2*a)
        r2 = (-b - (b**2-4*a*c)**0.5)/(2*a)
        return r1, r2
    else:
        print('a is zero, input again')
#%%
def add_end(L=[]):
    L.append('END')
    return L

#%%
#def person(name, age, **kws):
#    print('name:',name,'age:', age, 'others:',kws)
#
#extra={'city':'Beijing','job':'Engineer'}
#person('Jack', 24, city='Beijing',job='Engineer')
#person('Jack',24,**extra)
# function args
def f2(a, b, c=0, *, d, **kw):
    s = 'a=%s,b=%s,c=%s,d=%s,'
    t = []
    s1 = []
    for keys,vals in kw.items():
        s1.append(keys+'=%s')
        t.append(vals)
#        print(s1)
#        print(t)
    p = s+ ','.join(s1)
#    print(p)    
    print(p % (a,b,c,d,*t))
f2(1,4,d=6,e=7,f=8,g=20,h=22)    

#%%
# Hannor tower
def move(n, a, b, c):
    if n == 1:
        print('move', a, '-->', c)
        return
    move(n-1, a, c, b)
    print('move', a, '-->', c)
    move(n-1, b, a, c)

move(3, 'A', 'B', 'C')
#%%
with open('C:\\Python_Garage\\hello.txt','w') as f:
    f.write('hello python!')
with open('C:\\Python_Garage\\hello.txt') as f:    
    s=f.read(2)
    print(s)
#%%
import os
os.name    
import muser_rdraw as mr
import rdin as ri
from astropy.time import Time
#%%

pathname = r'C:\MUSER_DATA\data'
trange = ri.set_tr2db(pathname)
sfiles = ri.sel_file(pathname,trange)
ttick = 5
with open(sfiles[2],'rb') as f:
    gt1 = mr.gps_rd_m1(f,0)
    gt2 = mr.gps_rd_m1(f, 19199)
print(gt1,gt2)
dt = Time(gt2)-Time(gt1)
dtsec = dt.sec
print(dtsec)
#%%
def qround(num,step):
    a = step/2
    k = int(num/step)
    s = (k+1)*step if num-k*step >=a else k*step
    return s

#%%
from datetime import datetime
def log(message,when=None):
    when = datetime.now() if when==None else when
    print('%s %s' % (when, message))
#%%
from decimal import Decimal
rate = Decimal('1.45')
seconds = Decimal('222') # 3*60 + 42
cost = rate * seconds / Decimal('60')
print(cost)