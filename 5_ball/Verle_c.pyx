import Pair as p
import scipy.constants
import copy as cop

def new_cord (i):
    i.coord=i.coord+i.speed*1+0.5*i.a

def new_a(i, parts):
    G=scipy.constants.G
    i.last_a=cop.copy(i.a)
    i.a=p.Pair(0, 0)
    for j in parts:
        if abs(i.coord-j.coord)>=16:
            i.a=i.a+G*j.m*(j.coord-i.coord)/(abs(j.coord-i.coord)**3)
    i.a=i.a*10000000000

def new_vel (i):
    i.speed=i.speed+0.5*(i.a+i.last_a)*1

def update (parts):
    for i in parts:
        new_cord(i)
        new_a(i, parts)
        new_vel(i)
        i.time=i.time-1