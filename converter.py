import numpy as np

def read_poscar(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Read the header
    system_name = lines[0].strip()

    # Read the scaling factor
    scale = float(lines[1].strip())

    # Read the lattice vectors
    lattice_vectors = []
    for i in range(2, 5):
        lattice_vectors.append([float(x) for x in lines[i].split()])
    lattice_vectors = np.array(lattice_vectors) * scale

    # Read the element symbols and counts
    elements = lines[5].split()
    counts = [int(x) for x in lines[6].split()]

    # Read the coordinate style
    coord_style = lines[7].strip().lower()

    # Read the atomic coordinates
    coords = []
    symbols = []
    start = 8
    for i, count in enumerate(counts):
        end = start + count
        symbols.extend([elements[i]] * count)
        for j in range(start, end):
            coords.append([float(x) for x in lines[j].split()[:3]])
        start = end

    # If the coordinate style is fractional, convert to Cartesian
    if coord_style == 'direct' or coord_style == 'cartesian':
        pass
    elif coord_style == 'fractional':
        coords = np.dot(coords, lattice_vectors)
    else:
        raise ValueError('Unsupported coordinate style: {}'.format(coord_style))

    return system_name, lattice_vectors, symbols, coords

def write_cel(filename, system_name, lattice_vectors, symbols, coords):
    with open(filename, 'w') as f:
        f.write('Comment string\n')
        f.write('0 {:.4f} {:.4f} {:.4f} 90.0000 90.0000 90.0000\n'.format(
            lattice_vectors[0][0]/10, lattice_vectors[1][1]/10, lattice_vectors[2][2]/10))
        for symbol, coord in zip(symbols, coords):
            f.write('{} {:.6f} {:.6f} {:.6f} 1.000000 0.005100 0.000000 0.000000 0.000000\n'.format(
                symbol, coord[0], coord[1], coord[2]))
        f.write('*\n')


if __name__ == "__main__":
    poscar_file = 'BTO_[110]9[-11-4]_POSCAR'
    cel_file = 'BTO_[110]9[-11-4].cel'

    system_name, lattice_vectors, symbols, coords = read_poscar(poscar_file)
    write_cel(cel_file, system_name, lattice_vectors, symbols, coords)