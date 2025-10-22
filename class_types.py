class Node:
    def __init__ (self, id, x, y):
        self.id = id
        self.x = x
        self.y = y

class Element:
    def __init__ (self, id, nodes_id, jakobian_matrix):
        self.id = id
        self.nodes_id = nodes_id
        self.jakobian_matrix = jakobian_matrix
        
class Grid:
    def __init__ (self, nodes, elements):
        self.nodes = nodes
        self.elements = elements

class Jakobian:
    def __init__(self, J, J1, detJ):
        self.J = J
        self.J1 = J1
        self.detJ = detJ
