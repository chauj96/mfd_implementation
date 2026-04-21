import os
import numpy as np
from scipy.io import loadmat

# Mesh loader for different test cases (currently supports twoFaults)
# TO DO: add SPE11b and distorted unit cube case
 
def load_mesh(mesh_name):

    if mesh_name == "twoFaults":
        return load_twofault()

    else:
        raise ValueError(f"Unknown mesh: {mesh_name}")


def load_twofault():

    filepath = os.path.join("meshes", "twoFaults", "fault_mesh.mat")

    data = loadmat(filepath, simplify_cells=True)

    cell_struct = data["cell_struct"]
    face_struct = data["face_struct"]
    vertices = data["V3"]

    # Index conversion (from matlab to python)
    for c in cell_struct:
        c["faces"] = np.array(c["faces"]).astype(int) - 1

    # Size of the domain
    Lx = 1.0
    Ly = 1.0
    Lz = 0.4

    return cell_struct, face_struct, vertices, Lx, Ly, Lz