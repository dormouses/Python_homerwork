#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import Emitter as e
import Pair as p
import copy as cop
import scipy.constants
import thread 
import numpy as np
from scipy import integrate
from math import sqrt
import Verle_c

class Verle:
    def __init__(self):
        self.parts = []
        self.em = e.Emitter()
        self.dt = 1
        self.type = 0
        self.time = 0
        
    def make_part(self, x, y, vx, vy, m, t, col):
        self.parts.append(self.em.MakePart(x, y, vx, vy, m, t, col))
        
    def new_cord (self, i):
        i.coord=i.coord+i.speed*self.dt+0.5*i.a
    
    def new_a(self, i):
        G=scipy.constants.G
        i.last_a=cop.copy(i.a)
        i.a=p.Pair(0, 0)
        for j in self.parts:
            if abs(i.coord-j.coord)>=16:
                i.a=i.a+G*j.m*(j.coord-i.coord)/(abs(j.coord-i.coord)**3)
        i.a=i.a*10000000000
            
    def new_vel (self, i):
        i.speed=i.speed+0.5*(i.a+i.last_a)*self.dt
    
    def update_normal(self, i):
            self.new_cord(i)
            self.new_a(i)
            self.new_vel(i)
            i.time=i.time-1
    
    def kill(self):
        for i in self.parts:
            if i.time==0:
                self.parts.remove(i)
            
    def update_odeint(self, max_t, max_step):
        self.time=self.time %  max_step
        if self.time == 0:
            self.ode = self.count_odeint(max_t, max_step)  
        count = len(self.parts)
        for i in range(count):
            print '---'
            print self.time
            print 2*len(self.parts)+2*i
            print '---'
            self.parts[i].coord.x = self.ode[self.time][2*i]
            print '||||'
            print self.time
            print 2*len(self.parts)+2*i
            print '||||'
            self.parts[i].coord.y = self.ode[self.time][2*i+1]
            print '+++'
            print self.time
            print 2*len(self.parts)+2*i
            print '+++'
            self.parts[i].speed.x = self.ode[self.time][2*len(self.parts)+2*i]
            print '==='
            print self.time
            print 2*len(self.parts)+2*i
            print '==='
            self.parts[i].speed.y = self.ode[self.time][2*len(self.parts)+2*i+1]
        self.time=self.time+1;
     
    def count_odeint(self, max_t, max_step):
        G=scipy.constants.G
        t = np.linspace(0, max_t, max_step)
        count = len(self.parts)
        
        def Solve(y,t):
            result = np.zeros(4*count)
            for i in range(count):
                result[2*i] = y[2*count+2*i]
                result[2*i+1] = y[2*count+2*i+1]
                for j in range(count):
                    if j != i:
                        r = sqrt((y[2*j]-y[2*i])**2+(y[2*j+1]-y[2*i+1])**2)
                        if r >= 4:
                            result[2*count+2*i]+=G*10000000000*self.parts[j].m*(y[2*j]-y[2*i])/(r**3)
                            result[2*count+2*i+1]+=G*10000000000*self.parts[j].m*(y[2*j+1]-y[2*i+1])/(r**3)
            return result
            
        init = np.zeros(4*count)
        for i in range(count):
            init[2*i] = self.parts[i].coord.x
            init[2*i+1] = self.parts[i].coord.y
            init[2*count + 2*i] = self.parts[i].speed.x
            init[2*count + 2*i+1] = self.parts[i].speed.y
        return integrate.odeint(Solve, init, t)
   
           
    def update(self):
        self.kill()
        if self.type==0:
            for i in self.parts:
                self.update_normal(i)
        elif self.type==1:
            for i in self.parts:
                thread.start_new_thread(self.update_normal, (i,) )
        elif self.type==2:
            self.update_odeint(2, 2)
        elif self.type==3:
            Verle_c.update(self.parts)
            
                

        
