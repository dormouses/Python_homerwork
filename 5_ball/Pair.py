from math import sqrt

class Pair :
    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    def __add__(self, other):
        return Pair(self.x + other.x, self.y + other.y)
    
    def __sub__(self, other):
        return Pair(self.x - other.x, self.y - other.y)
    
    def __mul__(self, other):
        self.x = self.x * other
        self.y = self.y * other
        return self
    
    def __rmul__(self, other):
        self.x = self.x * other
        self.y = self.y * other
        return self
    
    def __div__(self, other):
        self.x = self.x / other
        self.y = self.y / other
        return self
    
    def __abs__(self):
        x=self.x
        y=self.y
        return sqrt(x*x + y*y)
        
    def __eq__ (self, other):
        return self.x==other.x and self.y==other.y