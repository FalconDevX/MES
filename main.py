from file_parse import GlobalData
from gauss_integrations import DerivativeCoordinates
from class_types import SystemOfEquations

if __name__ == "__main__":
    data = GlobalData("Test1_4_4.txt")
    print("N:",data.N)
    jacobian = DerivativeCoordinates(data.grid, data.Conductivity, data.N, data.BC, data.Alfa, data.Tot, data.Density, data.SpecificHeat)
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


    