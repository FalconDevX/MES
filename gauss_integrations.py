import numpy as np
from itertools import product
import inspect
from class_types import Node, Element, Grid, Jakobian
from utils import clean_near_zero

class GaussTable:
    def __init__(self, N : int):
        self.N = N
        self.nodes, self.weights = self.generate_gauss_table(N)

    def generate_gauss_table(self, N : int):
        if N == 1:
            x_k = [0.0]
            A_k = [2.0]
        elif N == 2:
            x_k = [-1/np.sqrt(3), 1/np.sqrt(3)]
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

#całkwoanie metoda gasussa
class GaussIntegral:
    def __init__(self, N: int, function):
        self.N = N 
        self.function = function
        self.dim = len(inspect.signature(function).parameters)
        self.gauss_table = GaussTable(self.N)

    def integrate(self):    
        weights = self.gauss_table.weights
        nodes = self.gauss_table.nodes
        integration_value = 0.0

        for points in product(range(self.N), repeat=self.dim):
            # print(points)
            w = 1.0
            coord = []
            for p in points:    
                w *= weights[p]
                coord.append(nodes[p])
            
            integration_value += w * self.function(*coord)
        return integration_value
    
#tabele pochodnych funkcji ksztaltu po ksi i eta
class DerivativeTable:
    def __init__(self):
        self.derivatives_ksi, self.derivatives_eta = self.generate_derivative_table()

    def  generate_derivative_table(self):
        derivatives_ksi = np.zeros((4, 4))
        derivatives_eta = np.zeros((4, 4))
        gauss_nodes = GaussTable(2).nodes

        gauss_points = list(product(gauss_nodes, repeat=2))
        
        print(gauss_points)

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

        for i, (ksi, eta) in enumerate(gauss_points):
            dNksi = dN_dksi(eta)
            dNeta = dN_deta(ksi)
            derivatives_ksi[i, :] = dNksi
            derivatives_eta[i, :] = dNeta

        # for i in range(derivatives_ksi.size):
        #     value = dN_dksi(gauss_points[i//4][1])
        #     derivatives_ksi.flat[i] = value[i%4]

        # for i in range(derivatives_eta.size):
        #     value = dN_deta(gauss_points[i//4][0])
        #     derivatives_eta.flat[i] = value[i%4]
            
        for i in derivatives_ksi:
            print(i)
            
        return derivatives_ksi, derivatives_eta

# obliczenia jakobianu 
class DerivativeCoordinates:
    def __init__ (self, grid : Grid, conductivity):
        der_table = DerivativeTable()
        self.conductivity = conductivity
        self.nodes = grid.nodes
        self.gauss_weights = GaussTable(2).weights
        self.der_table_eta = der_table.derivatives_eta
        self.der_table_ksi = der_table.derivatives_ksi
        self.elements = self.calculateJacobianMatrix(grid.elements, grid.nodes, self.der_table_eta, self.der_table_ksi, self.gauss_weights, self.conductivity)
        
    #obliczanie macierzy jacobiego, odwrotnej i wyznacznika
    def calculateJacobianMatrix(self, elements, nodes, der_table_eta, der_table_ksi, gauss_weights, conductivity):
    # for element in elements:
        element = elements[0]
        nodes_map = {n.id: n for n in nodes}
        x_coords = []
        y_coords = []

        element.jakobian = []
                
        for node_id in element.nodes_id:
            node = nodes_map[node_id]
            x_coords.append(node.x)
            y_coords.append(node.y)        
        H_local = []

        #liczenie 4 maceirzy jakobiego dla jednego elementu 
        for i in range(4):
            jakobian = Jakobian([],[],[])
            
            dx_dksi = clean_near_zero(sum(x_coords[j] * der_table_ksi[i][j] for j in range(4)))
            dy_dksi = clean_near_zero(sum(y_coords[j] * der_table_ksi[i][j] for j in range(4)))
            dx_deta = clean_near_zero(sum(x_coords[j] * der_table_eta[i][j] for j in range(4)))
            dy_deta = clean_near_zero(sum(y_coords[j] * der_table_eta[i][j] for j in range(4)))

            jakobian.J = np.array([[dx_dksi, dx_deta],[dy_dksi, dy_deta]])

            jakobian.detJ = np.linalg.det(jakobian.J)
            jakobian.J1 = np.linalg.inv(jakobian.J)
            element.jakobian.append(jakobian)
            der_eta_ksi_row_num = i
            H_local.append(self.calculateHMatrix(der_table_eta, der_table_ksi, jakobian, der_eta_ksi_row_num, conductivity))

        element.H = sum(H_local[i] * gauss_weights[0]*gauss_weights[1] for i in range(4))
        # return elements

    def calculateHMatrix(self, der_table_eta, der_table_ksi, jakobian: Jakobian, der_eta_ksi_row_num, conductivity):
        print("jakobian.J1", jakobian.J1)
        print("jakobian.detJ", jakobian.detJ)

        J1 = jakobian.J1  
        print("J1", J1)
        derivatives_x = np.zeros((4,1))
        derivatives_y = np.zeros((4,1))

        for i in range(4):
            dN_dksi = der_table_ksi[der_eta_ksi_row_num][i]
            dN_deta = der_table_eta[der_eta_ksi_row_num][i]
            dn_dx = J1[0][0] * dN_dksi + J1[0][1] * dN_deta
            derivatives_x[i] = dn_dx
            dn_dy = J1[1][0] * dN_dksi + J1[1][1] * dN_deta
            derivatives_y[i] = dn_dy

        print("Derivatives in x:")
        print(derivatives_x)
        print("Derivatives in y:")
        print(derivatives_y)

        dN_dx = derivatives_x.reshape(-1, 1)   
        dN_dy = derivatives_y.reshape(-1, 1)   

        dN_dx_T = dN_dx.T   
        dN_dy_T = dN_dy.T  

        #illoczyn wektora i wekotra transpozycji 
        dN_dx_prod = dN_dx @ dN_dx_T   
        dN_dy_prod = dN_dy @ dN_dy_T   
        #liczenie macierzy lokalnej dla jednego pc
        H_local = conductivity * (dN_dx_prod + dN_dy_prod) * jakobian.detJ
        print("H lokalne: ", H_local)
        return H_local


    # def print_jakobian(self):
    #     for element in self.elements:
    #         print(f"Element ID: {element.id}")
    #         for i, jakobian in enumerate(element.jakobian):
    #             print(f"  Gauss point {i+1}:")
    #             print(f"    J:\n{jakobian.J}")
    #             print(f"    detJ: {jakobian.detJ}")
    #             print(f"    J1:\n{jakobian.J1}")
    #             print()