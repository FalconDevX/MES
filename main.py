from file_parse import GlobalData
from gauss_integrations import DerivativeCoordinates
from class_types import SystemOfEquations

import numpy as np

def compute_t1(H, C, dtau, t0, P):
    H = np.array(H)
    C = np.array(C)

    print("H: ", H)
    print("C: ", C)
    print("dtau: ", dtau)
    print("t0: ", t0)
    print("P: ", P)

    # Zrób z t0 wektor jeśli podano skalara
    if np.isscalar(t0):
        t0 = np.full(H.shape[0], t0)

    t0 = np.array(t0).reshape(-1)
    P  = np.array(P).reshape(-1)

    A = H + C / dtau
    B = (C / dtau) @ t0 - P

    
    print("A: ", A)
    print("B: ", B)
    print("t1: ", np.linalg.solve(A, B))

    return np.linalg.solve(A, B)

if __name__ == "__main__":
    data = GlobalData("Test1_4_4.txt")
    print("N:",data.N)
    jacobian = DerivativeCoordinates(data.grid, data.Conductivity, data.N, data.BC, data.Alfa, data.Tot, data.Density, data.SpecificHeat)
    
    t1 = compute_t1(jacobian.H_global, jacobian.C_global, data.SimulationStepTime, data.InitialTemp, jacobian.P_global)
    # jacobian.print_jakobian()

    # system_of_equations = SystemOfEquations(jacobian.H_global, jacobian.P_global)
    # system_of_equations.solve()
    # print("Temperature: ", system_of_equations.t)


    # for element in data.elements:
    #     print(f"Element ID: {element.id}, Nodes IDs: {element.nodes_id}")

    # print("Wszystkie węzły")
    # for node in data.nodes:
    #     print(f"ID={node.id}, x={node.x}, y={node.y}")

    # print("\nElementy i ich węzły")
    # for element in data.elements:
    #     print(f"Element {element.id}: {element.nodes_id}")

    # print("\nID węzłów poszczególnych elementów")
    # for element in data.elements:
    #     print(f"Element E{element.id} ma węzły o ID: {element.nodes_id}")


    