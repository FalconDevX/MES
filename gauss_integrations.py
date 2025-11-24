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

#funkcje kształtu dla macierzy Hbc(całak powierzchniowa)
class ShapeFunctions:
    def N_functions(self, ksi, eta):
        return np.array([
        0.25 * (1 - ksi) * (1 - eta),
        0.25 * (1 + ksi) * (1 - eta),
        0.25 * (1 + ksi) * (1 + eta),
        0.25 * (1 - ksi) * (1 + eta),
    ])

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

    def generate_derivative_table(self):
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
        # print("ksi: ")
        # for i in derivatives_ksi:
        #     print(i)
        # #eta
        # print("eta: ")
        # print()
        # for i in derivatives_eta:
        #     print(i)
        return derivatives_ksi, derivatives_eta

# obliczenia jakobianu 
class DerivativeCoordinates:
    def __init__ (self, grid : Grid, conductivity, N, BC, alfa):
        der_table = DerivativeTable(N)
        self.alfa = alfa
        self.N = N
        self.BC = BC
        self.conductivity = conductivity
        self.nodes = grid.nodes
        self.gauss_weights = GaussTable(N).weights
        self.der_table_eta = der_table.derivatives_eta
        self.der_table_ksi = der_table.derivatives_ksi
        self.elements = self.calculateJacobianMatrix(grid.elements, grid.nodes, self.der_table_eta, self.der_table_ksi, self.gauss_weights, self.conductivity, self.BC, self.alfa)
        self.H_global = self.agregateHmatrix(self.elements, self.nodes)
        
        np.set_printoptions(precision=4, suppress=True, linewidth=200)
        print(self.H_global)

    #obliczanie macierzy jacobiego, odwrotnej i wyznacznika
    def calculateJacobianMatrix(self, elements, nodes, der_table_eta, der_table_ksi, gauss_weights, conductivity, BC, alfa):

        for element in elements:
            #mapping node id to node coordinates
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

            # print("elementy brzegowe")

            ##Obkliczanie macierzy Hbc

            #tworzenie punktów lokalnych brzegowych 
            
            gauss_nodes = GaussTable(self.N).nodes
            gauss_weights = GaussTable(self.N).weights

            #gauss_points = [(ksi, eta) for eta in gauss_table for ksi in gauss_table]

            ksi_eta_edge_points = {
                # dla pierwszej krawedzi[( -1/√3 , -1 ),(  1/√3 , -1 )]
                1: [(ksi, -1) for ksi in gauss_nodes],   
                2: [(1,  eta) for eta in gauss_nodes],   
                3: [(ksi,  1) for ksi in gauss_nodes],   
                4: [(-1, eta) for eta in gauss_nodes]    
            }

            # for i in ksi_eta_edge_points:
            #     print("ksi_eta_edge_points: ", ksi_eta_edge_points[i])

            #liczenie Hbc w jednej petli z H
            element.Hbc = np.zeros((4,4))
            #sprawdzanie czy element jest elementem brzegowym
            edges_indexes = [(0,1), (1,2), (2,3), (3,0)]
            element_nodes = element.nodes_id
            
            #id1, id2, id_edge
            boundary_edges = []
            # 0 0 1
            # 1 1 2
            # 2 2 3
            # 3 3 0
            for id, (i,j) in enumerate(edges_indexes):
                #jeśli węzły są w BC
                if element_nodes[i] in BC and element_nodes[j] in BC:
                    #id + 1 bo id nie zaczyna się od 0
                    boundary_edges.append((element_nodes[i], element_nodes[j], id + 1))

            shape_functions = ShapeFunctions()

            #część właściwa obliczania macierzy Hbc
            for bd_edge in boundary_edges:
                # print("bd_edge: ", bd_edge)
                #uzyskiwanie punktów ksi i eta na podstawie numeru krawędzi
                edge_id = bd_edge[2]
                edge_points = ksi_eta_edge_points[edge_id]       
                # print("edge_points: ", edge_points)
                #pc1 , pc2 dwa punkty na tej samej krawędzi 
                pc1 = np.array(shape_functions.N_functions(edge_points[0][0], edge_points[0][1])).reshape(4,1)
                pc2 = np.array(shape_functions.N_functions(edge_points[1][0], edge_points[1][1])).reshape(4,1)
                # print("pc1: ", pc1)
                # print("pc2: ", pc2)
                #iloczyn wektorowy
                H_pc1 = pc1 @ pc1.T   
                H_pc2 = pc2 @ pc2.T   

                #liczenie jakobianu na krawedzi 
                id1, id2, edge_id = bd_edge
                #uzyskiwanie współrzędnych węzłów na podstawie numeru węzła
                p1 = nodes_map[id1]
                p2 = nodes_map[id2]

                dx = p2.x - p1.x
                dy = p2.y - p1.y
                J_edge = np.sqrt(dx*dx + dy*dy) / 2

                # print("H_pc1: ", H_pc1 * gauss_weights[0] * alfa * J_edge)
                # print("H_pc2: ", H_pc2 * gauss_weights[1] * alfa * J_edge)

                Hbc_local = alfa * (H_pc1 * gauss_weights[0] + H_pc2 * gauss_weights[1]) * J_edge
                # print("Hbc_local: ", Hbc_local)
                element.Hbc += Hbc_local
            element.H += element.Hbc
            # print("--------------------------------")
            # print("element.id: ", element.id)
            # print("element.Hbc: ", element.Hbc)
        return elements

    def calculateHMatrix(self, der_table_eta, der_table_ksi, jakobian: Jakobian, der_eta_ksi_row_num, conductivity, element: Element):
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

    def agregateHmatrix(self, elements, nodes):
        all_nodes_num = len(nodes)
        H_global = np.zeros((all_nodes_num, all_nodes_num))

        for element in elements:
            ids = element.nodes_id
            H_local = element.H  

            for i in range(4):
                for j in range(4):
                    i_global = ids[i] - 1
                    j_global = ids[j] - 1
                    H_global[i_global, j_global] += H_local[i,j]
        return H_global

    def print_jakobian(self):
        GREEN = "\033[32m"      
        BLUE = "\033[34m"
        CYAN = "\033[36m"
        MAGENTA = "\033[35m"
        RESET = "\033[0m"

        for element in self.elements:
            print(f"{GREEN} Element ID:{RESET} {element.id}")
            print(f"    {GREEN}H{RESET}:\n{element.H}")
            print(f"    {GREEN}Hbc{RESET}:\n{element.Hbc}")
            # for i, jakobian in enumerate(element.jakobian):
            #     print(f"  {CYAN}Gauss point {i+1}:{RESET}")
            #     print(f"    {GREEN}J{RESET}:\n{jakobian.J}")
            #     print(f"    {BLUE}detJ{RESET}: {jakobian.detJ}")
            #     print(f"    {MAGENTA}J1{RESET}:\n{jakobian.J1}")
            #     print()
            for i, H in enumerate(element.H_local):
                print(f"  {CYAN}H local{i+1}:{RESET}\n{H}")
            # for i, der_table in enumerate(element.der_table):
            #     print(f"  {CYAN}Derivatives table {i+1}:{RESET}\n{der_table}")
