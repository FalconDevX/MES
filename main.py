from file_parse import GlobalData
from mes import MesMatrixElements
from class_types import SystemOfEquations

import numpy as np

test1 = "Test1_4_4.txt"
test2 = "Test2_4_4_MixGrid.txt"
test3 = "Test3_31_31_kwadrat.txt"

def compute_t1(H, C, dtau, t0, P, return_matrices=False):
    H = np.array(H)
    C = np.array(C)
    
    #t0 jako wektor kolumnowy
    if np.isscalar(t0):
        t0 = np.full((H.shape[0], 1), t0)
    else:
        t0 = np.array(t0).reshape(-1, 1)
    
    P = np.array(P).reshape(-1, 1)
    
    # Macierz A = H + C/dtau
    A = H + C / dtau
    
    B = (C / dtau) @ t0 + P
    
    # Rozwiązanie układu równań
    t1 = np.linalg.solve(A, B)
    
    if return_matrices:
        return t1, A, B
    return t1

def simulation(H_global, C_global, P_global, T0, SimulationTime, dt):
    T_old = T0.copy()
    num_steps = int(SimulationTime / dt)
    
    # Ustaw opcje wyświetlania NumPy, aby pokazywać wszystkie wiersze
    np.set_printoptions(threshold=np.inf, linewidth=np.inf)
    
    print("\nMax and min temperature in each step")
    print(f"{'Time[s]':<12} {'MinTemp[C]':<15} {'MaxTemp[C]':<15}")
    print("-" * 42)
    
    for step in range(num_steps):
        
        T_new, A, B = compute_t1(H_global, C_global, dt, T_old, P_global, return_matrices=True)
        
        #Macierz A
        # print(f"\nInteration {step}")
        # print(f"H Matrix ([H]+[C]/dT)")
        # print(A)
        
        #Macierz B
        # print(f"\nP_Vector (({{P}}+{{[C]/dT}}*{{TO}}))")
        # print(B.flatten())
        T_new = compute_t1(H_global, C_global, dt, T_old, P_global)
        
        #Obliczanie min i max temperatury
        min_temp = np.min(T_new)
        max_temp = np.max(T_new)
        current_time = (step + 1) * dt
        
        print(f"\n{current_time:<12.0f} {min_temp:<15.6f} {max_temp:<15.6f}")
        
        T_old = T_new
    
    return T_new

if __name__ == "__main__":
    data = GlobalData(test3)
    print("N:",data.N)
    mes = MesMatrixElements(data.grid, data.Conductivity, data.N, data.BC, data.Alfa, data.Tot, data.Density, data.SpecificHeat)
    
    #mes.print_jakobian()
    # zamiana temperatury początkowej na wektor globalny
    num_nodes = len(data.nodes)
    T0 = np.full((num_nodes, 1), data.InitialTemp)
    
    simulation(mes.H_global, mes.C_global, mes.P_global, T0, data.SimulationTime, data.SimulationStepTime)
    # mes.print_jakobian()

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


    