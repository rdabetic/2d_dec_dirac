import matplotlib.pyplot as plt
import matplotlib.tri as tri
import meshio
import sys

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script_name.py <fname (without .msh)>")
        sys.exit(1)

    fname = sys.argv[1]

    mesh = meshio.read(fname + ".vtk")
    T = tri.Triangulation(mesh.points[:, 0],
                          mesh.points[:, 1],
                          mesh.cells[0].data)
    plt.triplot(T)
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(fname + ".pdf", dpi = 300)

