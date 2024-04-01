def read_xyz(filename):
    atoms = []
    coordinates = []
    with open(filename) as xyz:
        n_atoms = int(xyz.readline().strip())
        xyz.readline()  # skip the blank line

        for _ in range(2, n_atoms + 2):
            line = xyz.readline()
            atom, x, y, z = line.split()
            atoms.append(atom)
            coordinates.append([float(x), float(y), float(z)])
    return atoms, coordinates
