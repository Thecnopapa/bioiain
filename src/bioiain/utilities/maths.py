import numpy as np
from types import GeneratorType
from ..biopython.base import BiopythonOverlayClass

def find_com(atoms:list|GeneratorType|np.ndarray|BiopythonOverlayClass) -> list[float]:
    """
    Find the center of mass of a list of atoms. All atoms weight the same.
    :param atoms: List (or generator) of Atom objects. Also accepts as atoms lists and np.arrays, where the first 3
    items are the coordinates.
    :return: List of coordinates of the center of mass.
    """
    x = 0
    y = 0
    z = 0

    if isinstance(atoms, BiopythonOverlayClass):
        atoms = atoms.get_atoms()

    if isinstance(atoms, GeneratorType):
        atoms = list(atoms)

    for atom in atoms:
        if isinstance(atom, np.ndarray) or type(atom) == list:
            x += atom[0]
            y += atom[1]
            z += atom[2]
        else:
            x += atom.coord[0]
            y += atom.coord[1]
            z += atom.coord[2]
    x /= len(atoms)
    y /= len(atoms)
    z /= len(atoms)
    return [x, y, z]



def vector(b, e):
    if len(b) == 3:
        x, y, z = b
        X, Y, Z = e
        return (X - x, Y - y, Z - z)
    elif len(b) == 2:
        x, y = b
        X, Y = e
        return (X - x, Y - y)
    return None


def d2(p0, p1, root=False):
    if root:
        return np.sqrt((p0[0] - p1[0]) ** 2 + (p0[1] - p1[1]) ** 2 + (p0[2] - p1[2]) ** 2)
    else:
        return (p0[0] - p1[0]) ** 2 + (p0[1] - p1[1]) ** 2 + (p0[2] - p1[2]) ** 2


def length(v):
    if len(v) == 3:
        x, y, z = v
        return np.sqrt(x * x + y * y + z * z)
    elif len(v) == 2:
        x, y = v
        return np.sqrt(x * x + y * y)
    return None


def distance(p0, p1):
    return d2(p0, p1, root=True)