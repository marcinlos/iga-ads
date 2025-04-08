import numpy as np


def data_reader(lines: list, K: int) -> np.ndarray:
    '''
    This function reads the porosity/permeability data from the Poros.data/Perm.data file and returns a numpy array of the porosity/permeability data for the given plane K.

    Parameters
    ----------
    lines : list
        List of strings containing the lines of the Poros.data/Perm.data file.
    K : int
        The plane number for which the data is to be read. The original source files had K ranging from 1 to 3.

    Returns
    -------
    data_plane: np.ndarray
        A numpy array of the porosity/permeability data for the given plane K.

    Raises
    ------
    None
    '''

    data_plane = np.zeros((25, 25))
    K_start = False

    for line in lines:
        if "Plane K = {}".format(K) in line:
            K_start = True
            continue
        if "Plane K = {}".format(K+1) in line:
            break
        if K_start:
            if line.strip().startswith("I ="):
                temp_map = line.strip().split()[2:] # this will skip everything up to I = and get the I indices
                a = int(temp_map[0]) - 1
                b = int(temp_map[-1])
            if line.strip().startswith("J="):
                row_indx = int(line.strip().split()[1]) - 1 # this will get the J index
                row_data = line.strip().split()[2:] # Skip the J= part and get the data
                row_data = [float(i) for i in row_data]
                data_plane[row_indx, a:b] = row_data

    return data_plane.T # Transpose the data to have x-y coordinates


# read porosity data from 'Poros.data'

with open('Poros.data') as f:
    lines = f.readlines()

poros_k1 = data_reader(lines, 1)
poros_k2 = data_reader(lines, 2)
poros_k3 = data_reader(lines, 3)

# read permeability data from 'Perm.data'

with open('Perm.data') as f:
    lines = f.readlines()

perm_k1 = data_reader(lines, 1)
perm_k2 = data_reader(lines, 2)
perm_k3 = data_reader(lines, 3)



def data_ascii_wrapper(data: np.ndarray, filename: str, a_x: float = None, b_x: float = None, a_y: float = None, b_y: float = None) -> None:
    '''
    This function writes the data to an ascii file in the format: {x}\t{y}\t{poros/perm value}.

    Parameters:
    ----------
    data: np.ndarray
        The porosity or permeability data to be written to the file.
        Must be of the same format as the data_reader() function output.
        I.e. 2D numpy array with the first dimension being the x-axis and the second dimension being the y-axis and the values being the porosity/permeability values.
    filename: str
        The name of the file to write the data to.
    a_x: float
        The starting value of the x-axis. If None, the starting value will be 0.
    b_x: float
        The ending value of the x-axis. If None, the ending value will be the length of the x-axis.
    a_y: float
        The starting value of the y-axis. If None, the starting value will be 0.
    b_y: float
        The ending value of the y-axis. If None, the ending value will be the length of the y-axis.

    Returns:
    -------
    None
    '''

    a_x = a_x or 0
    b_x = b_x or data.shape[0]
    step_x = (b_x - a_x) / data.shape[0]

    a_y = a_y or 0
    b_y = b_y or data.shape[1]
    step_y = (b_y - a_y) / data.shape[1]

    with open(filename, 'w') as f:
        f.write(f"{data.shape[0]}\t{data.shape[1]}\n")
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                f.write(f"{a_x + i * step_x:.2f}\t{a_y + j * step_y:.2f}\t{data[i, j]:.3f}\n")
        f.close()


# wrap the porosity and permeability data to ascii files

data_ascii_wrapper(poros_k1, 'porosity_k1.data')
data_ascii_wrapper(poros_k2, 'porosity_k2.data')
data_ascii_wrapper(poros_k3, 'porosity_k3.data')

data_ascii_wrapper(perm_k1, 'permeability_k1.data')
data_ascii_wrapper(perm_k2, 'permeability_k2.data')
data_ascii_wrapper(perm_k3, 'permeability_k3.data')

# interpolate the data to a finer grid

from scipy.interpolate import RegularGridInterpolator

def data_interpolator(data: np.ndarray, a_x: float = 25, b_x: float = 25, nsteps_x: int = 100, a_y: float = 25, b_y: float = 25, nsteps_y: int = 100) -> np.ndarray:
    '''
    This function interpolates the data to a finer grid using RegularGridInterpolator from scipy.

    Parameters:
    ----------
    data: np.ndarray
        The porosity or permeability data to be interpolated.
        Must be of the same format as the data_reader() function output.
        I.e. 2D numpy array with the first dimension being the x-axis and the second dimension being the y-axis and the values being the porosity/permeability values.
    a_x: float
        The starting value of the x-axis. Default is 25.
    b_x: float
        The ending value of the x-axis. Default is 25.
    nsteps_x: int
        The number of steps to interpolate the x-axis. Default is 100.
    a_y: float
        The starting value of the y-axis. Default is 25.
    b_y: float
        The ending value of the y-axis. Default is 25.
    nsteps_y: int
        The number of steps to interpolate the y-axis. Default is 100.

    Returns:
    -------
    interp_data: np.ndarray
        The interpolated data on the finer grid.
    '''

    x = np.linspace(0, data.shape[0], data.shape[0])
    y = np.linspace(0, data.shape[1], data.shape[1])

    interpolator = RegularGridInterpolator((x, y), data, method='linear')

    x_new = np.linspace(0, a_x, nsteps_x)
    y_new = np.linspace(0, b_y, nsteps_y)
    xv, yv = np.meshgrid(x_new, y_new)
    points = np.array([xv.flatten(), yv.flatten()]).T
    interp_data = interpolator(points).reshape(nsteps_x, nsteps_y)

    return interp_data.T

# interpolate the porosity and permeability data to a finer grid

# poros_k1_interp = data_interpolator(poros_k1)
# poros_k2_interp = data_interpolator(poros_k2)
# poros_k3_interp = data_interpolator(poros_k3)

# perm_k1_interp = data_interpolator(perm_k1)
# perm_k2_interp = data_interpolator(perm_k2)
# perm_k3_interp = data_interpolator(perm_k3)

# wrap the interpolated data to ascii files

# data_ascii_wrapper(poros_k1_interp, 'porosity_k1_interp.data', a_x=0, b_x=25, a_y=0, b_y=25)
# data_ascii_wrapper(poros_k2_interp, 'porosity_k2_interp.data', a_x=0, b_x=25, a_y=0, b_y=25)
# data_ascii_wrapper(poros_k3_interp, 'porosity_k3_interp.data', a_x=0, b_x=25, a_y=0, b_y=25)

# data_ascii_wrapper(perm_k1_interp, 'permeability_k1_interp.data', a_x=0, b_x=25, a_y=0, b_y=25)
# data_ascii_wrapper(perm_k2_interp, 'permeability_k2_interp.data', a_x=0, b_x=25, a_y=0, b_y=25)
# data_ascii_wrapper(perm_k3_interp, 'permeability_k3_interp.data', a_x=0, b_x=25, a_y=0, b_y=25)
