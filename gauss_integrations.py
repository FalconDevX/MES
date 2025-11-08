import numpy as np
from itertools import product
import inspect
from class_types import Element, Grid, Jakobian

class GaussTable:
    def __init__(self, N : int):
        self.N = N
        self.nodes, self.weights = self.generate_gauss_table(N)

    def generate_gauss_table(self, N : int):
        """
            Returns Gauss–Legendre points and weights for the given N (1–4).
        """
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

#całkwowanie metoda gasussa
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
    
    def __init__(self, N):
        self.N = N
        self.derivatives_ksi, self.derivatives_eta = self.generate_derivative_table()

    def  generate_derivative_table(self):
        """
            Calculate derivatives table from N shape functions and local vairables (ksi, eta) for every integral point
            arg: N - number of integral points scheme
            returns: array of integrated shape functions after local poinst (ksi, eta)
        """
        gauss_nodes = GaussTable(self.N).nodes

        gauss_points = [(ksi, eta) for eta in gauss_nodes for ksi in gauss_nodes]

        derivatives_ksi = np.zeros((len(gauss_points), 4))
        derivatives_eta = np.zeros((len(gauss_points), 4))

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
            # dNksi = np.round(dN_dksi(eta), 6)
            # dNeta = np.round(dN_deta(ksi), 6)
            derivatives_ksi[i, :] = dN_dksi(eta)
            derivatives_eta[i, :] = dN_deta(ksi)

        # for i in range(derivatives_ksi.size):
        #     value = dN_dksi(gauss_points[i//4][1])
        #     derivatives_ksi.flat[i] = value[i%4]

        # for i in range(derivatives_eta.size):
        #     value = dN_deta(gauss_points[i//4][0])
        #     derivatives_eta.flat[i] = value[i%4]
        
        #ksi
        print("ksi: ")
        for i in derivatives_ksi:
            print(i)
        #eta
        print("eta: ")
        print()
        for i in derivatives_eta:
            print(i)
        return derivatives_ksi, derivatives_eta

# obliczenia jakobianu 
class DerivativeCoordinates:
    def __init__ (self, grid : Grid, conductivity, N):
        der_table = DerivativeTable(N)
        self.N = N
        self.conductivity = conductivity
        self.nodes = grid.nodes
        self.gauss_weights = GaussTable(N).weights
        self.der_table_eta = der_table.derivatives_eta
        self.der_table_ksi = der_table.derivatives_ksi
        self.elements = self.calculateJacobianMatrix(grid.elements, grid.nodes, self.der_table_eta, self.der_table_ksi, self.gauss_weights, self.conductivity)
        
    #obliczanie macierzy jacobiego, odwrotnej i wyznacznika
    def calculateJacobianMatrix(self, elements, nodes, der_table_eta, der_table_ksi, gauss_weights, conductivity):
        """
            Calculate jacobian matrix, transposed and determinant for every element in one iteration
            arg: class contructor
            returns: array of elements
        """
        for element in elements:
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
            for i in range(len(der_table_ksi)):
                jakobian = Jakobian([],[],[])
                
                dx_dksi = sum(x_coords[j] * der_table_ksi[i][j] for j in range(4))
                dy_dksi = sum(y_coords[j] * der_table_ksi[i][j] for j in range(4))
                dx_deta = sum(x_coords[j] * der_table_eta[i][j] for j in range(4))
                dy_deta = sum(y_coords[j] * der_table_eta[i][j] for j in range(4))

                jakobian.J = np.array([[dx_dksi, dy_dksi],[dx_deta, dy_deta]])

                jakobian.detJ = np.linalg.det(jakobian.J)
                jakobian.J1 = np.linalg.inv(jakobian.J)
                element.jakobian.append(jakobian)
                der_eta_ksi_row_num = i
                H_local.append(self.calculateHMatrix(der_table_eta, der_table_ksi, jakobian, der_eta_ksi_row_num, conductivity, element))
            element.H_local = H_local
            #iloczyn kartezjański dla różnych wag dla 2d
            gauss_weights_2d = list(product(gauss_weights, repeat=2)) 
            
            #para elementów bo ksi i eta
            element.H = sum(H_local[i] * w_pair[0] * w_pair[1] for i, w_pair in enumerate(gauss_weights_2d))
        return elements

    def calculateHMatrix(self, der_table_eta, der_table_ksi, jakobian: Jakobian, der_eta_ksi_row_num, conductivity, element: Element):
        """
            Calculate H local matrix for every node
            arg: derivatives tables for eta and ksi, jacobian object, row number from 
        """

        J1 = jakobian.J1
        derivatives_x = np.zeros((4, 1))
        derivatives_y = np.zeros((4, 1))

        #4 bo dla 2D tylko ksi i eta: 4 x dN_dksi, 4 x dN_deta, 
        for i in range(4):
            #wiersz zgodny z "i" z petli wywołującej calculateHmatrix czyli liczba wierszy
            dN_dksi = der_table_ksi[der_eta_ksi_row_num][i]
            dN_deta = der_table_eta[der_eta_ksi_row_num][i]
            #wzór na elementy pochodznych loklanych po funkcji kształtu
            dn_dx = J1[0][0] * dN_dksi + J1[0][1] * dN_deta
            derivatives_x[i] = dn_dx
            dn_dy = J1[1][0] * dN_dksi + J1[1][1] * dN_deta
            derivatives_y[i] = dn_dy
            element.der_table.append((derivatives_x.copy(), derivatives_y.copy()))

        # print("Derivatives in x:")
        # print(derivatives_x)
        # print("Derivatives in y:")
        # print(derivatives_y)

        dN_dx = derivatives_x.flatten()
        dN_dy = derivatives_y.flatten()

        H_local = np.zeros((4,4))
        for i in range(4):
            for j in range(4):
                #pełny wzór na mceirz lokalna
                H_local[i,j] = (dN_dx[i]*dN_dx[j] + dN_dy[i]*dN_dy[j]) * conductivity * abs(jakobian.detJ)
        return H_local


    def print_jakobian(self):
        GREEN = "\033[32m"      
        BLUE = "\033[34m"
        CYAN = "\033[36m"
        MAGENTA = "\033[35m"
        RESET = "\033[0m"

        for element in self.elements:
            print(f"{GREEN} Element ID:{RESET} {element.id}")
            print(f"    {GREEN}H{RESET}:\n{element.H}")
            for i, jakobian in enumerate(element.jakobian):
                print(f"  {CYAN}Gauss point {i+1}:{RESET}")
                print(f"    {GREEN}J{RESET}:\n{jakobian.J}")
                print(f"    {BLUE}detJ{RESET}: {jakobian.detJ}")
                print(f"    {MAGENTA}J1{RESET}:\n{jakobian.J1}")
                print()
            for i, H in enumerate(element.H_local):
                print(f"  {CYAN}H local{i+1}:{RESET}\n{H}")
            # for i, der_table in enumerate(element.der_table):
            #     print(f"  {CYAN}Derivatives table {i+1}:{RESET}\n{der_table}")
