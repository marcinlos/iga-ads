import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

# import matplotlib animation
from matplotlib.animation import FuncAnimation

import argparse
parser = argparse.ArgumentParser(description='Create a movie of CO2 sequestration simulation')
parser.add_argument('--n_iter', type=int, default=100, help='Number of iterations')
parser.add_argument('--t_step', type=float, default=1, help='Time step, s')
parser.add_argument('--fps', type=int, default=3, help='Frames per second')
args = parser.parse_args()

n_iter = args.n_iter
t_step = args.t_step 
fps = args.fps

plt.rcParams['font.size'] = 16
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['figure.titlesize'] = 20

# load data
data_dict = {}
data_dict['p'] = {}
data_dict['s'] = {}
for i in tqdm(range(n_iter)):
    data_dict['p'][i] = np.loadtxt('p.out_' + str(i * 10) + '.data')
    data_dict['s'][i] = np.loadtxt('s.out_' + str(i * 10) + '.data')

def heatmap_data(data_dict, label, i):
    data = data_dict[label][i]

    x = data[:, 0]
    y = data[:, 1]
    vals = data[:, 2]

    x_unique = np.unique(x)
    y_unique = np.unique(y)
    heatmap_data = np.zeros((len(x_unique), len(y_unique)))

    for i in range(len(x)):
        x_idx = np.where(x_unique == x[i])[0][0]
        y_idx = np.where(y_unique == y[i])[0][0]
        heatmap_data[x_idx, y_idx] = vals[i]

    return heatmap_data.T

fig, ax = plt.subplots(1, 2, figsize=(36, 12))

p_heatmap = ax[0].imshow(heatmap_data(data_dict, 'p', 0), cmap='viridis', aspect='equal', animated=True)
s_heatmap = ax[1].imshow(heatmap_data(data_dict, 's', 0), cmap='viridis', aspect='equal', animated=True)
p_colorbar = fig.colorbar(p_heatmap, ax=ax[0])
s_colorbar = fig.colorbar(s_heatmap, ax=ax[1])

ax[0].invert_yaxis()
ax[1].invert_yaxis()

def animate(i):
    p_heatmap.set_array(heatmap_data(data_dict, 'p', i))
    s_heatmap.set_array(heatmap_data(data_dict, 's', i))

    vmin_p = np.min(heatmap_data(data_dict, 'p', i))
    vmax_p = np.max(heatmap_data(data_dict, 'p', i))
    vmin_s = np.min(heatmap_data(data_dict, 's', i))
    vmax_s = np.max(heatmap_data(data_dict, 's', i))

    p_heatmap.set_clim(vmin=vmin_p, vmax=vmax_p)
    s_heatmap.set_clim(vmin=vmin_s, vmax=vmax_s)
    
    p_colorbar.update_normal(p_heatmap)
    s_colorbar.update_normal(s_heatmap)

    fig.suptitle('Time: ' + str(i * t_step) + ' s')

    return p_heatmap, s_heatmap


ax[0].set_xlabel('X-axis')
ax[0].set_ylabel('Y-axis')
ax[0].set_title('Pressure')
ax[1].set_xlabel('X-axis')
ax[1].set_ylabel('Y-axis')
ax[1].set_title('Saturation')

ani = FuncAnimation(fig, animate, frames=tqdm(range(n_iter)), interval=100)
ani.save('heatmap.mp4', writer='ffmpeg', fps=fps)


