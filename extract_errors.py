import re
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 16})
import sys
import numpy as np

def extract_errors(filename):
    l2_norm_errors = []
    sobolev_norm_errors = []
    mesh_widths = []
    
    # Regular expressions to match the required fields
    l2_norm_pattern = re.compile(r'Error in DEC L2 norm:\s*([\d.e+-]+)')
    sobolev_norm_pattern = re.compile(r'Error in DEC Sobolev norm:\s*([\d.e+-]+)')
    mesh_width_pattern = re.compile(r'Mesh-width:\s*([\d.e+-]+)')
    
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            # Search for L2 norm error
            l2_match = l2_norm_pattern.search(line)
            if l2_match:
                l2_norm_errors.append(float(l2_match.group(1)))
            # Search for Sobolev norm error
            sobolev_match = sobolev_norm_pattern.search(line)
            if sobolev_match:
                sobolev_norm_errors.append(float(sobolev_match.group(1)))
            # Search for Mesh-width
            mesh_width_match = mesh_width_pattern.search(line)
            if mesh_width_match:
                mesh_widths.append(float(mesh_width_match.group(1)))

    return np.array(l2_norm_errors), \
            np.array(sobolev_norm_errors), \
            np.array(mesh_widths)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script_name.py <load_fname (without .txt)> \
              <optional: reference order>")
        sys.exit(1)

    load_fname = sys.argv[1]
    order = 1
    if len(sys.argv) >= 3:
        order = int(sys.argv[2])
    
    l2_errors, sobolev_errors, mesh_widths = \
        extract_errors(load_fname + ".txt")

    plt.loglog(mesh_widths, l2_errors, label = r"DEC $L^2$", marker = '+')
    plt.loglog(mesh_widths, sobolev_errors, label = r"DEC $H\Lambda$", marker = '+')
    plt.loglog(mesh_widths, 
               (mesh_widths * sobolev_errors[0] **  (1 / order) /
                mesh_widths[0]) ** order, 
               label = r"$\mathcal{O}$" + f"$(h^{{{order}}})$",
               ls = ":")
    plt.xlabel("Mesh-Width")
    plt.ylabel("Error norms")
    plt.legend(loc = 'upper left')
    plt.tight_layout()
    plt.savefig(load_fname + ".pdf", dpi = 300)
