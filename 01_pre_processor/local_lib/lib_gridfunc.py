import sys
import numpy as np

def uniform(x, y, z, line):
    """Distribute grid points uniformly."""
    if line['dir'] in ['X', 'x']:
        xyz = x
    elif line['dir'] in ['Y', 'y']:
        xyz = y
    elif line['dir'] in ['Z', 'z']:
        xyz = z
    else:
        raise ValueError(f"[Error] Invalid direction '{line['dir']}' for uniform distribution.")

    delta = (line['coord_f'] - line['coord_i']) / (line['index_f'] - line['index_i'])

    for i in range(line['index_i'], line['index_f'] + 1):
        xyz[i] = line['coord_i'] + (i - line['index_i']) * delta

    print(f"{line['dir'].lower()}-dir., Uniform, {line['coord_i']:.3f} ~ {line['coord_f']:.3f}, "
          f"Interval is {delta:.6f}.")

def geometric(x, y, z, line):
    """Distribute grid points using geometric progression."""
    if line['dir'] in ['X', 'x']:
        xyz = x
    elif line['dir'] in ['Y', 'y']:
        xyz = y
    elif line['dir'] in ['Z', 'z']:
        xyz = z
    else:
        raise ValueError(f"[Error] Invalid direction '{line['dir']}' for geometric distribution.")

    factor1 = line['factor1']  # Expansion (1) or Compression (0)
    factor2 = line['factor2']  # Geometric Ratio (> 1)

    if factor2 <= 1:
        raise ValueError("[Error] For geometric progression, ratio must be larger than 1.")
    if factor1 != 1:
        factor2 = 1 / factor2

    init_delta = (line['coord_f'] - line['coord_i']) * (1 - factor2) / (1 - factor2 ** (line['index_f'] - line['index_i']))
    final_delta = init_delta * factor2 ** (line['index_f'] - line['index_i'] - 1)

    for i in range(line['index_i'], line['index_f'] + 1):
        xyz[i] = line['coord_i'] + init_delta * (1 - factor2 ** (i - line['index_i'])) / (1 - factor2)

    print(f"{line['dir'].lower()}-dir., Geometric, {line['coord_i']:.3f} ~ {line['coord_f']:.3f}, "
          f"Interval from {init_delta:.6f} to {final_delta:.6f}.")

def hypertan(x, y, z, line):
    """Distribute grid points using hyperbolic tangent progression."""
    if line['dir'] in ['X', 'x']:
        xyz = x
    elif line['dir'] in ['Y', 'y']:
        xyz = y
    elif line['dir'] in ['Z', 'z']:
        xyz = z
    else:
        raise ValueError(f"[Error] Invalid direction '{line['dir']}' for hyperbolic tangent distribution.")

    factor1 = line['factor1']  # Expansion (1), Compression (0), or Symmetric (0.5)
    factor2 = line['factor2']  # Gamma value

    if factor1 != 0:
        for i in range(line['index_i'], line['index_f'] + 1):
            prop = (i - line['index_i']) / (line['index_f'] - line['index_i'])
            xyz[i] = line['coord_i'] + (line['coord_f'] - line['coord_i']) * factor1 * \
                     (1 - np.tanh(factor2 * (factor1 - prop)) / np.tanh(factor2 * factor1))
    else:
        factor1 = 1
        for i in range(line['index_i'], line['index_f'] + 1):
            prop = (i - line['index_i']) / (line['index_f'] - line['index_i'])
            xyz[i] = line['coord_i'] + (line['coord_f'] - line['coord_i']) * factor1 * \
                     (1 - np.tanh(factor2 * prop) / np.tanh(factor2 * factor1))

    init_delta = xyz[line['index_i'] + 1] - xyz[line['index_i']]
    mid_delta = xyz[int((line['index_i'] + line['index_f']) / 2) + 1] - xyz[int((line['index_i'] + line['index_f']) / 2)]
    final_delta = xyz[line['index_f']] - xyz[line['index_f'] - 1]

    print(f"{line['dir'].lower()}-dir., Hypertan, {line['coord_i']:.3f} ~ {line['coord_f']:.3f}, "
          f"Interval from {init_delta:.6f} through {mid_delta:.6f} to {final_delta:.6f}.")
