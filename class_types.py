class Node:
    def __init__ (self, id, x, y):
        self.id = id
        self.x = x
        self.y = y

class Element:
    def __init__ (self, id, nodes_id):
        self.id = id
        self.nodes_id = nodes_id

class Grid:
    def __init__ (self, nodes, elements):
        self.nodes = nodes
        self.elements = elements

