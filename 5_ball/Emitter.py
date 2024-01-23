import random as rand
import Pair as pair
import Part as p

class Emitter:

    def MakePart(self, x, y, vx, vy, m, time, col):
        coord = pair.Pair(x, y)
        speed = pair.Pair(vx, vy)
        return p.Part(coord, speed, m, time, col)