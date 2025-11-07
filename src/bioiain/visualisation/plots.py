import sys, os


from ..utilities.logging import log

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

mpl_colours = 'blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan'

try:
    mpl.use('QtAgg')
except:
    try:
        mpl.use('TkAgg')
    except:
        mpl.use('Agg')










def fig3D(entity=None, preset=None, fig=None, ax=None):
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


    return fig, ax