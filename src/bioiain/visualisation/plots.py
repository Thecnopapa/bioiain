import sys, os, math
from .. import TEMP_FOLDER

from ..utilities.logging import log

import numpy as np


import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d



# if os.environ.get("MPLBACKEND", None) is not None:
#     try:
#         mpl.use('QtAgg')
#     except:
#         try:
#             mpl.use('TkAgg')
#         except:
#             mpl.use('Agg')

log(1, f"MPL backend: {mpl.get_backend()}")

mpl_colours = ('blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan')
pymol_colours = ('green', 'cyan', 'red', 'yellow', 'violet','blue',
               'salmon', 'lime', 'pink', 'slate', 'magenta', 'orange', 'marine',
               'olive', 'purple', 'teal', 'forest', 'firebrick', 'chocolate',
               'wheat', 'white', 'grey')

def plasma(value, scale=256, as_hex=False, as_pymol_hex=False, alpha:float|None=None):
    cm = mpl.colormaps["plasma"]
    if scale <= 0:
        scale = 1
    #print(value, scale)
    value = round((float(value) / float(scale)) * 256)
    #print(value)
    col = cm(value)
    col = process_rgba(col, alpha=alpha, as_hex=as_hex, as_pymol_hex=as_pymol_hex)
    return col

class Plasma():
    def __init__(scale=256, as_hex=False, as_pymol_hex=False, alpha:float|None=None):
        self.scale=scale
        self.as_hex=as_hex
        self.as_pymol_hex=as_pymol_hex
        self.alpha=alpha

    def __call__(self, scale=None, as_hex=None, as_pymol_hex=None, alpha:float|None=None):
        if scale is None:
            scale=self.scale
        if as_hex is None:
            as_hex=self.as_hex
        if as_pymol_hex is None:
            as_pymol_hex



def traffic(value, scale=256, as_hex=False, as_pymol_hex=False, alpha:float|None=None):
    cm = mpl.colormaps["RdYlGn"]
    if scale <= 0:
        scale = 1
    #print(value, scale)
    value = round((float(value) / float(scale)) * 256)
    #print(value)
    col = cm(value)
    col = process_rgba(col, alpha=alpha, as_hex=as_hex, as_pymol_hex=as_pymol_hex)
    return col

def process_rgba(col, alpha=None, as_hex=False, as_pymol_hex=False):
    if alpha is not None:
        col = list(col)
        col[3] = alpha
    if as_hex or as_pymol_hex:
        from matplotlib.colors import rgb2hex
        col = rgb2hex(col)
    if as_pymol_hex:
        col = col.replace("#", "0x")
    return col





def close(fig):
    plt.close(fig)

def show(**kwargs):
    plt.show(**kwargs)


def grid2D(rows, columns, height=5, width=5, as_grid=False):
    log(2, f"Creating {rows}x{columns} (rxc) grid...")
    fig, grid_axes = plt.subplots(rows, columns, figsize=(columns*width, rows*height))
    if as_grid:
        return fig, grid_axes
    axes = []
    try:
        for row in grid_axes:
            axes.extend(row)
    except TypeError:
        axes = grid_axes

    return fig, axes

def fig2D(**kwargs):
    fig = plt.figure(**kwargs)
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    return fig, ax


def line(start, end):
    return list(zip(start, end))



def fig3D(entity:any=None, preset:str=None,
          fig:mpl.pyplot.Figure=None,
          ax:mpl.pyplot.Axes=None) -> list[mpl.pyplot.Figure|mpl.pyplot.Axes]:
    """
    Initialises a 3-Dimensional Matplotlib plot, optionally using a preset.
    :param entity: (Optional) The entity to plot, might be required in some presets.
    :param preset: Name of preset to use.
    :param fig: Use this figure instead of new one.
    :param ax: Use these axes instead of new ones.
    :return: Figures and Axes generated.
    """
    if fig is None:
        fig = plt.figure()
    if ax is None:
        ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')

    if preset == "crystal-frac":
        ax.set_xticks([-1, 0, 1, 2])
        ax.set_yticks([-1, 0, 1, 2])
        ax.set_zticks([-1, 0, 1, 2])

        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")

        ax.axes.set_xlim(-1,2)
        ax.axes.set_ylim(-1,2)
        ax.axes.set_zlim(-1,2)

    else:
        #log("warning", "3D Plot preset ({}) not found".format(preset))
        pass


    return [fig, ax]





class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        arrow_prop_dict = dict(mutation_scale=20, arrowstyle='-|>', color='k', shrinkA=0, shrinkB=0)
        kwargs = arrow_prop_dict | kwargs
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))

        return np.min(zs)



def plot_heatmap(matrix, show=False, filename=None):
    fig = plt.figure()
    ax = fig.add_subplot()

    ax.imshow(matrix)
    if filename is not None:
        plt.savefig(filename)
    if show:
        plt.show(block=True)


def mpl3D_to_gif(
    fig, 
    axes, 
    name = "animation.gif", 
    folder=None, 
    save_path=None,
    total_frames=360, total_d = 360, duration = 15
    ):


    import io, PIL
    from ..utilities.utilities_old import ProgressBar

    if type(axes) not in (list, tuple):
        axes = [axes]

    if folder is None:
        folder = TEMP_FOLDER

    if save_path is not None:
        name = os.path.basename(save_path)
        folder = os.path.dirname(save_path)

    if not name.endswith(".gif"):
        name += ".gif"
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, name)


    f_duration = duration * 1000 / total_frames

    dpf = total_d / total_frames

    log(3, "Animating: {} Duration: {}s, degrees/frames: {}/{}".format(name, duration, total_d, total_frames))


    progress = ProgressBar(total_frames, silent=True)
    images = []
    for frame in range(total_frames):
        for ax in axes:
            ax.view_init(azim=frame*dpf)
        buf = io.BytesIO()
        fig.savefig(buf)
        buf.seek(0)

        images.append(PIL.Image.open(buf))
        progress.add(info="{}/{}".format(frame, total_frames))
    images[0].save(
        path,
        append_images=images[1:],
        duration=f_duration,  # duration of each frame in milliseconds
        loop=1,  # loop forever
        save_all=True,
    )
    log(4,"Saving to:", path)
    return path

def voxels3d(value_grid, count_grid=None, show_plot=False, title=None, shrink=False, gif_path=None, ax=None, **kwargs):

    def explode_cube(data):
        size = np.array(data.shape)*2
        data_e = np.zeros(size - 1, dtype=data.dtype)
        data_e[::2, ::2, ::2] = data
        return data_e

    print()
    log(3, f"Plotting voxels... ({title})")

    value_grid = np.copy(value_grid)

    if ax is None:
        fig, ax = fig3D()
    else:
        fig = None
    np.set_printoptions(threshold=sys.maxsize)
    # print(value_grid)
    max_val = abs(value_grid.reshape(value_grid.shape[-1] ** 3).max() - value_grid.reshape(value_grid.shape[-1] ** 3).min())
    norm_grid = value_grid - value_grid.reshape(value_grid.shape[-1] ** 3).min()

    log(4, "Scale:", max_val)
    if title is None:
        title = "voxels"
    ax.set_title(f"{title} s={max_val:.2f}")

    # cube = np.indices([size, size, size])

    if count_grid is None:
        log(4, "Plotting entire cube (no count_grid)")
        count_grid = value_grid.reshape(value_grid.shape)
        count_grid.fill(1)
    else:
        count_grid = np.copy(count_grid)
    cube = (count_grid > 0) & (count_grid > 0) & (count_grid > 0)
    # print(cube)


    color_vector = np.vectorize(plasma)
    colors = np.array(color_vector(norm_grid, scale=max_val, as_hex=True))
    # print(colors)

    log(4, "Cube:", cube.shape)
    log(4, "Colors:", colors.shape)



    if shrink:
        #filled = np.ones(cube.shape)
        filled = explode_cube(cube)
        colors = explode_cube(colors)
        #ecolors_2 = explode(edgecolors)
        x, y, z = np.indices(np.array(filled.shape) + 1).astype(float) // 2
        x[0::2, :, :] += 0.05
        y[:, 0::2, :] += 0.05
        z[:, :, 0::2] += 0.05
        x[1::2, :, :] += 0.95
        y[:, 1::2, :] += 0.95
        z[:, :, 1::2] += 0.95

        ax.voxels(x, y, z, filled, facecolors=colors, alpha=0.5)

    else:
        ax.voxels(cube, facecolors=colors, alpha=0.5)

    # ax.set_box_aspect((size, size, size))
    ax.set_aspect('equal')

    if fig is not None:
        if gif_path is not None:
            mpl3D_to_gif(fig, ax, save_path=gif_path)

        if show_plot:
            show()
        return fig, ax
    else:
        return ax
