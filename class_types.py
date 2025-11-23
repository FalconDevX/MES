class Node:
    def __init__ (self, id, x, y, BC):
        self.id = id
        self.x = x
        self.y = y
        self.BC = BC

class Element:
    def __init__ (self, id, nodes_id, jakobian, H, H_local, der_table, Hbc):
        self.id = id
        self.nodes_id = nodes_id
        self.jakobian = jakobian
        self.H = H
        self.H_local = H_local
        self.der_table = der_table
        self.Hbc = Hbc
        
class Grid:
    def __init__ (self, nodes, elements):
        self.nodes = nodes
        self.elements = elements

class Jakobian:
    def __init__(self, J, J1, detJ):
        self.J = J
        self.J1 = J1
        self.detJ = detJ

class Surface:
    def __init__(self, N):
        self.N = N