from class_types import Node, Element, Grid

class GlobalData:
    def __init__(self, filename):
        self.filename = filename
        self.N = None
        self.SimulationTime = None
        self.SimulationStepTime = None
        self.Conductivity = None
        self.Alfa = None
        self.Tot = None
        self.InitialTemp = None
        self.Density = None
        self.SpecificHeat = None

        self.nodes = []
        self.elements = []
        self.BC = [] 

        self.read_file()
        self.grid = Grid(self.nodes, self.elements)

    def read_file(self):
        with open(self.filename, 'r') as file:
            lines = file.readlines()

        mode = None
        for line in lines:
            line = line.strip()
            if not line:
                continue

            if line.startswith("IntegralScheme"):
                self.N = float(line.split()[1])
            elif line.startswith("SimulationTime"):
                self.SimulationTime = float(line.split()[1])
            elif line.startswith("SimulationStepTime"):
                self.SimulationStepTime = float(line.split()[1])
            elif line.startswith("Conductivity"):
                self.Conductivity = float(line.split()[1])
            elif line.startswith("Alfa"):
                self.Alfa = float(line.split()[1])
            elif line.startswith("Tot"):
                self.Tot = float(line.split()[1])
            elif line.startswith("InitialTemp"):
                self.InitialTemp = float(line.split()[1])
            elif line.startswith("Density"):
                self.Density = float(line.split()[1])
            elif line.startswith("SpecificHeat"):
                self.SpecificHeat = float(line.split()[1])
            elif line.startswith("*Node"):
                mode = "nodes"
                continue
            elif line.startswith("*Element"):
                mode = "elements"
                continue
            elif line.startswith("*BC"):
                mode = "bc"
                continue

            if mode == "nodes":
                parts = [x.strip() for x in line.split(',') if x.strip()]
                if len(parts) == 3:
                    ID, x, y = parts
                    self.nodes.append(Node(int(ID), float(x), float(y)))

            elif mode == "elements":
                parts = [x.strip() for x in line.split(',') if x.strip()]
                if len(parts) == 5:
                    ID, n1, n2, n3, n4 = parts
                    self.elements.append(Element(int(ID), [int(n1), int(n2), int(n3), int(n4)], [], [], [], []))

            elif mode == "bc":
                parts = [x.strip() for x in line.split(',') if x.strip()]
                self.BC.extend(map(int, parts))
