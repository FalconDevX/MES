# MES - 2D Transient Heat Conduction FEM Solver

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![NumPy](https://img.shields.io/badge/NumPy-1.20+-green.svg)](https://numpy.org/)
[![Matplotlib](https://img.shields.io/badge/Matplotlib-3.0+-orange.svg)](https://matplotlib.org/)
[![License](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)

A **Finite Element Method (FEM)** solver for 2D transient heat conduction with quadrilateral (4-node) elements, Gauss quadrature integration, and Robin (convection) boundary conditions.

---

## Features

- **2D heat equation** - Transient conduction with capacity matrix **C**, conductivity matrix **H**, and load vector **P**
- **Quadrilateral mesh** - 4-node isoparametric elements (DC2D4)
- **Gauss integration** - 1-4 integration points (configurable via input)
- **Boundary conditions** - Robin (convection) on edges: $q = \alpha(T - T_{ot})$
- **Implicit time stepping** - Backward Euler: $(H + C/\Delta\tau)\,T^{n+1} = C/\Delta\tau\,T^n + P$
- **Input format** - Simple text mesh files (nodes, elements, material/simulation parameters)
- **Visualization** - Optional heatmap of temperature field via Matplotlib

---

## Installation

### Requirements

- Python 3.8+
- NumPy
- Matplotlib (optional, for heatmaps)

### Setup

```bash
git clone https://github.com/your-username/MES.git
cd MES
pip install numpy matplotlib
```

Or with a requirements file:

```bash
pip install -r requirements.txt
```

Create `requirements.txt` if needed:

```
numpy>=1.20
matplotlib>=3.0
```

---

## Usage

Run the solver with a mesh file (default in code is `Test2_4_4_MixGrid.txt`):

```bash
python main.py
```

To use another mesh, change the test file in `main.py`:

```python
# In main.py
test1 = "Test1_4_4.txt"
test2 = "Test2_4_4_MixGrid.txt"
test3 = "Test3_31_31_kwadrat.txt"
test4 = "mesh_100x100.txt"

# Then:
data = GlobalData(test2)  # or test1, test3, test4
```

Output: global **H** and **C** matrices, and per time step the min/max temperature.

To enable heatmap plotting at each step, uncomment in `main.py`:

```python
plot_heatmap(data.grid, T_new)
```

---

## Input File Format

Mesh and simulation parameters are read from a text file with the following structure:

| Keyword           | Description                          | Example   |
|-------------------|--------------------------------------|-----------|
| `IntegralScheme`  | Gauss integration points (1–4)       | `4`       |
| `SimulationTime`  | Total simulation time [s]            | `500`     |
| `SimulationStepTime` | Time step [s]                     | `50`      |
| `Conductivity`    | Thermal conductivity [W/(m·K)]       | `25`      |
| `Alfa`            | Convection coefficient [W/(m²·K)]    | `300`     |
| `Tot`             | Ambient temperature [°C]             | `1200`    |
| `InitialTemp`     | Initial temperature [°C]             | `100`     |
| `Density`         | Density [kg/m³]                      | `7800`    |
| `SpecificHeat`    | Specific heat [J/(kg·K)]             | `700`     |

Then:

- **\*Node** - Lines: `id, x, y`
- **\*Element, type=DC2D4** - Lines: `id, n1, n2, n3, n4` (node IDs)
- **\*BC** - Comma-separated node IDs on the boundary (Robin BC)

Example:

```
IntegralScheme 4
SimulationTime 500
SimulationStepTime 50
Conductivity 25
Alfa 300
Tot 1200
InitialTemp 100
Density 7800
SpecificHeat 700
Nodes number 16
Elements number 9
*Node
  1, 0.1, 0.005
  2, 0.0666, 0.005
  ...
*Element, type=DC2D4
  1, 1, 2, 6, 5
  ...
*BC
1, 2, 3, 4, 5, 8, 9, 12, 13, 14, 15, 16
```

---

## Project Structure

```
MES/
├── main.py          # Entry point: load mesh, build matrices, run simulation
├── mes.py           # FEM: Gauss tables, shape functions, H/C/P assembly, BC
├── file_parse.py    # Parse mesh and global data from input file
├── class_types.py   # Node, Element, Grid, Jakobian, SystemOfEquations
├── heatmap.py       # Temperature heatmap visualization
├── utils.py         # Utilities
├── README.md
├── Test1_4_4.txt
├── Test2_4_4_MixGrid.txt
├── Test3_31_31_kwadrat.txt
└── mesh_100x100.txt
```

---

## Theory (short)

- **PDE:** $\rho c_p \frac{\partial T}{\partial t} = k \nabla^2 T + \dot{q}$, with convection on boundary.
- **Weak form** leads to $C \dot{T} + H T = P$.
- **Time discretization:** backward Euler → $(H + C/\Delta\tau) T^{n+1} = (C/\Delta\tau) T^n + P$.
- **H** and **Hbc** (boundary) from element conductivity and convection; **C** from capacity; **P** from boundary flux. All assembled from 4-node quadrilateral elements with Gauss quadrature.

---

## License

MIT License - see [LICENSE](LICENSE) for details.
