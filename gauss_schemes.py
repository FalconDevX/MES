import numpy as np
from itertools import product


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
    def __init__(self, N: int, function, dim: int):
        self.N = N
        self.function = function
        self.dim = dim
        self.gauss_table = GaussTable(N)

    def integrate(self):    
        weights = self.gauss_table.weights
        nodes = self.gauss_table.nodes
        integration_value = 0.0

        for points in product(range(self.N), repeat=self.dim):
            w = 1.0
            coord = []
            for p in points:
                w *= weights[p]
                coord.append(nodes[p])
            
            integration_value += w * self.function(*coord)
        return integration_value