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



def vector(b:list, e:list) -> list[float]|None:
    """
    Calculates a vector from given coordinates of 2 or 3 dimensions, otherwise returns None.
    :param b: First set of coordinates (start point).
    :param e: Second set of coordinates (end point).
    :return: Vector (e - b)
    """
    assert len(b) == len(e)
    if len(b) == 3:
        x, y, z = b
        X, Y, Z = e
        return [X - x, Y - y, Z - z]
    elif len(b) == 2:
        x, y = b
        X, Y = e
        return [X - x, Y - y]
    return None


def d2(p0:list, p1:list, root=False) -> float:
    """
    Calculates square distance between two points. Optionally return the square root of the square distance.
    :param p0: First set of coordinates.
    :param p1: Second set of coordinates.
    :param root: Whether to return the square root of the square distance.
    :return: Square distance (or root square distance).
    """
    if root:
        return np.sqrt((p0[0] - p1[0]) ** 2 + (p0[1] - p1[1]) ** 2 + (p0[2] - p1[2]) ** 2)
    else:
        return (p0[0] - p1[0]) ** 2 + (p0[1] - p1[1]) ** 2 + (p0[2] - p1[2]) ** 2


def length(v:list) -> float:
    """
    Calculates the length of a vector in space of 2 or 3 dimensions.
    :param v: Vector, can be generated using vector()
    :return: Length
    """
    if len(v) == 3:
        x, y, z = v
        return np.sqrt(x * x + y * y + z * z)
    elif len(v) == 2:
        x, y = v
        return np.sqrt(x * x + y * y)
    return None


def distance(p0, p1) -> float:
    """
    Returns the absolute distance between two points in space. Use d2() for square distance.
    :param p0: First set of coordinates.
    :param p1: Second set of coordinates.
    :return: Absolute distance.
    """
    return d2(p0, p1, root=True)



def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    if any(v): #if not all zeros then
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        return np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))

    else:
        return np.eye(3) #cross of all zeros only occurs on identical directions