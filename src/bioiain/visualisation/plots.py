import sys, os


from ..utilities.logging import log

import numpy as np


import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


from sklearn.decomposition import PCA
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sb

mpl_colours = 'blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan'

try:
    mpl.use('QtAgg')
except:
    try:
        mpl.use('TkAgg')
    except:
        mpl.use('Agg')










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
        log("warning", "3D Plot preset ({}) not found".format(preset))


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







# Confusion matrices
def plot_confusion(preds, labels, title, score=None, classes=None):
    try:
        log(1, "Plotting confusion...")
        cm = confusion_matrix(labels, preds)
        if classes is None:
            classes = list(set(labels))
        plt.figure(figsize=(1*len(classes) ,1*len(classes)))
        sb.heatmap(cm, annot=True, fmt="d", cmap="Blues", xticklabels=classes, yticklabels=classes)
        plt.xlabel("Predicted")
        plt.ylabel("True")
        if score is not None:
            title = f"{title}_S={score:.2f}"
        plt.title(title)
        os.makedirs("figs", exist_ok=True)
        path = f"figs/{title}.confusion.png"
        plt.savefig(path)
        plt.close()
        return cm, path
    except Exception as e:
        print(e)