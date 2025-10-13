from file_parse import GlobalData

if __name__ == "__main__":
    data = GlobalData("Test1_4_4.txt")

    print("Wszystkie węzły")
    for node in data.nodes:
        print(f"ID={node.id}, x={node.x:.8f}, y={node.y:.8f}")

    print("\nElementy i ich węzły")
    for element in data.elements:
        print(f"Element {element.id}: {element.nodes_id}")

    print("\nID węzłów poszczególnych elementów")
    for element in data.elements:
        print(f"Element E{element.id} ma węzły o ID: {element.nodes_id}")


    