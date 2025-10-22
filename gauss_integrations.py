import numpy as np
from itertools import product
import inspect

class GaussTable:
    def __init__(self, N : int):
        self.N = N
        self.nodes, self.weights = self.generate_gauss_table(N)

    def generate_gauss_table(self, N : int):
        if N == 1:
            x_k = [0.0]
            A_k = [2.0]
        elif N == 2:
            x_k = [1/np.sqrt(3), -1/np.sqrt(3)]
            A_k = [1.0, 1.0]
        elif N == 3:
            x_k = [np.sqrt(3/5), 0.0, -np.sqrt(3/5)]
            A_k = [5/9, 8/9, 5/9]
        elif N == 4:
            x_k = [np.sqrt(3/7 + 2/7 * np.sqrt(6/5)),
                    np.sqrt(3/7 - 2/7 * np.sqrt(6/5)),
                    -np.sqrt(3/7 - 2/7 * np.sqrt(6/5)),
                    -np.sqrt(3/7 + 2/7 * np.sqrt(6/5))]
            A_k = [(18 - np.sqrt(30)) / 36,
                    (18 + np.sqrt(30)) / 36,
                    (18 + np.sqrt(30)) / 36,
                    (18 - np.sqrt(30)) / 36]
        else:
            raise ValueError("Wrong argument N. Supported values are 1, 2, 3, 4.")

        return x_k, A_k
    
class GaussIntegral:
    def __init__(self, N: int, function):
        self.N = N + 1
        self.function = function
        self.dim = len(inspect.signature(function).parameters)
        self.gauss_table = GaussTable(self.N)

    def integrate(self):    
        weights = self.gauss_table.weights
        nodes = self.gauss_table.nodes
        integration_value = 0.0

        for points in product(range(self.N), repeat=self.dim):
            print(points)
            w = 1.0
            coord = []
            for p in points:    
                w *= weights[p]
                coord.append(nodes[p])
            
            integration_value += w * self.function(*coord)
        return integration_value
    
class DerivativeTable:
    def __init__(self):
        self.derivatives_dksi, self.derivatives_deta = self.generate_derivative_table()

    def  generate_derivative_table(self):
        derivatives_dksi = np.zeros((4, 4))
        derivatives_deta = np.zeros((4, 4))
        gauss_nodes = GaussTable(2).nodes

        gauss_points = list(product(gauss_nodes, repeat=2))

        def dN_dksi(eta):
            return [
                -0.25 * (1 - eta),
                0.25 * (1 - eta),
                0.25 * (1 + eta),
                -0.25 * (1 + eta)
            ]

        def dN_deta(ksi):
            return [
                -0.25 * (1 - ksi),
                -0.25 * (1 + ksi),
                0.25 * (1 + ksi),
                0.25 * (1 - ksi)
            ]

        for i in range(derivatives_dksi.size):
            value = dN_dksi(gauss_points[i//4][1])
            derivatives_dksi.flat[i] = value[i%4]

        for i in range(derivatives_deta.size):
            value = dN_deta(gauss_points[i//4][0])
            derivatives_deta.flat[i] = value[i%4]
            
        return derivatives_deta, derivatives_dksi
        
class DerivativeCoordinates:
    def __init__ (self, elements, nodes):
        self.elements = elements
        self.nodes = nodes
        self.der_table = DerivativeTable().derivatives_deta
        self.der_table_ksi = DerivativeTable().derivatives_dksi
        self.calculateJacobianMatrix(elements, nodes)

    def calculateJacobianMatrix(self, elements, nodes):
        for element in elements:
            nodes_map = {n.id: n for n in nodes}
            for node_id in element.nodes_id:
                node = nodes_map[node_id]
                print(f"Element ID: {element.id}, Node ID: {node.id}, x: {node.x}, y: {node.y}")