import numpy as np
import matplotlib.pyplot as plt

# I/O utilities for visualization and reporting:
# - export mesh and data to VTU (ParaView)
# - print solver error tables and plot the error

def write_vtu(filename, V3, cell_struct, face_struct, cellData, cellDataName):

    nCells = len(cell_struct)
    nPts = V3.shape[0]

    VTK_POLYHEDRON = 42

    connectivity = []
    offsets = []
    types = [VTK_POLYHEDRON] * nCells

    faces_all = []
    faceoffsets = []

    off_conn = 0
    off_face = 0

    for c in range(nCells):

        fids = np.array(cell_struct[c]["faces"]).astype(int)

        vids = []
        for f in fids:
            vids.extend(face_struct[f]["verts"])

        vids = np.unique(vids)
        vids0 = vids - 1  

        connectivity.extend(vids0)
        off_conn += len(vids0)
        offsets.append(off_conn)

        rec = [len(fids)]

        for f in fids:
            v = np.array(face_struct[f]["verts"]) - 1
            rec.extend([len(v)])
            rec.extend(v.tolist())

        faces_all.extend(rec)
        off_face += len(rec)
        faceoffsets.append(off_face)

    with open(filename, "w") as fid:

        fid.write('<?xml version="1.0"?>\n')
        fid.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        fid.write('<UnstructuredGrid>\n')
        fid.write(f'<Piece NumberOfPoints="{nPts}" NumberOfCells="{nCells}">\n')

        # Points
        fid.write('<Points>\n')
        fid.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        for p in V3:
            fid.write(f"{p[0]} {p[1]} {p[2]}\n")
        fid.write('</DataArray>\n</Points>\n')

        # Cells
        fid.write('<Cells>\n')

        fid.write('<DataArray type="Int32" Name="connectivity" format="ascii">\n')
        fid.write(" ".join(map(str, connectivity)))
        fid.write('\n</DataArray>\n')

        fid.write('<DataArray type="Int32" Name="offsets" format="ascii">\n')
        fid.write(" ".join(map(str, offsets)))
        fid.write('\n</DataArray>\n')

        fid.write('<DataArray type="UInt8" Name="types" format="ascii">\n')
        fid.write(" ".join(map(str, types)))
        fid.write('\n</DataArray>\n')

        fid.write('<DataArray type="Int32" Name="faces" format="ascii">\n')
        fid.write(" ".join(map(str, faces_all)))
        fid.write('\n</DataArray>\n')

        fid.write('<DataArray type="Int32" Name="faceoffsets" format="ascii">\n')
        fid.write(" ".join(map(str, faceoffsets)))
        fid.write('\n</DataArray>\n')

        fid.write('</Cells>\n')

        # Cell data
        fid.write(f'<CellData Scalars="{cellDataName}">\n')
        fid.write(f'<DataArray type="Int32" Name="{cellDataName}" format="ascii">\n')
        fid.write(" ".join(map(str, cellData)))
        fid.write('\n</DataArray>\n</CellData>\n')

        fid.write('</Piece>\n</UnstructuredGrid>\n</VTKFile>\n')

    # print(f"Wrote {filename}")

def print_pressure_err(results):

    print("\n=== Solver Results ===")
    print(f"{'tol':>10} | {'rel error':>12} | {'abs error':>12}")
    print("-"*40)

    for tol, rel_err, abs_err in results:
        print(f"{tol:10.1e} | {rel_err:12.3e} | {abs_err:12.3e}")


def plot_pressure_err(results):

    tol = np.array([r[0] for r in results])
    rel_err = np.array([r[1] for r in results])

    plt.figure()
    plt.loglog(tol, rel_err, '-o', label="Relative Error")
    plt.loglog(tol, tol, '--', label="y = tol")

    plt.xlabel("Tolerance")
    plt.ylabel("Relative Flux Error")
    plt.grid(True)
    plt.title("Flux Relative Error vs Tolerance")

    plt.legend(loc="upper left", frameon=True, fontsize=12)

    plt.show()