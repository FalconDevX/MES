from file_parse import GlobalData
from gauss_integrations import DerivativeCoordinates

if __name__ == "__main__":
    data = GlobalData("Test2_4_4_MixGrid.txt")
    jacobian = DerivativeCoordinates(data.elements, data.nodes)
    jacobian.print_jakobian()

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


    