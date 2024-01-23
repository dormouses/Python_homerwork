import Pair as p

class Part:
    def __init__(self, pos, speed, m, time, col):
        self.coord = pos
        self.speed = speed
        self.m = m
        self.time = time
        self.a= p.Pair(0,0)
        self.last_a = p.Pair(0,0)
        self.col=col
        self.r=m//50
        
    def __eq__ (self, other):
        return self.coord==other.coord