from .space_groups import dictio_space_groups
import Bio.PDB as bp
import numpy as np
from ..utilities import vector


def get_operation(key:int,op_n:int) -> dict:
    """
    Fetch translation and rotation parameters for a given symmetry operation in a given space group.
    :param key: Space group key
    :param op_n: Operation number
    :return: Translation / Rotation dict
    """
    return dictio_space_groups[key]["symops"][op_n]


def coord_add(coord:list[float], deltas:list[float], subtract:bool = False) -> list:
    """
    Add or subtract distances to a coordinate.
    :param coord: Coordinate (X, Y, Z)
    :param deltas: Distances (X, Y, Z)
    :param subtract: True to subtract
    :return: New coordinate
    """
    x, y, z = coord
    if subtract:
        nx = x - deltas[0]
        ny = y - deltas[1]
        nz = z - deltas[2]
    else:
        nx = x + deltas[0]
        ny = y + deltas[1]
        nz = z + deltas[2]
    return [nx, ny, nz]


def coord_operation(coord:list[float], key:int, op_n:int, distance:list[float] = (0,0,0)) -> list:
    """
    Perform a symmetry operation on a coordinate.
    :param coord: Coordinate (X, Y, Z)
    :param key: Space group key
    :param op_n: Operation number
    :param distance: (optional) Distance to add (X, Y, Z)
    :return: New coordinate
    """
    rotation = dictio_space_groups[key]
    operation = rotation["symops"][op_n]
    rot = operation["rot"]
    tra = operation["tra"]

    x, y, z = coord

    nx = (rot[0][0] * x) + (rot[0][1] * y) + (rot[0][2] * z) + tra[0]+distance[0]
    ny = (rot[1][0] * x) + (rot[1][1] * y) + (rot[1][2] * z) + tra[1]+distance[1]
    nz = (rot[2][0] * x) + (rot[2][1] * y) + (rot[2][2] * z) + tra[2]+distance[2]

    return [nx, ny, nz]


def coord_operation_entity(entity:bp.Entity.Entity, key:int, op_n:int, distance:list[float] = (0,0,0)) -> bp.Entity.Entity:
    """
    Perform a symmetry operation on an Entity in-place.
    :param entity: Entity to perform symmetry operation on
    :param key: Space group key
    :param op_n: Operation number
    :param distance: (optional) Distance to add (X, Y, Z)
    :return: Original Entity transformed
    """
    for atom in entity.get_atoms():
        if atom.is_disordered() > 0:
            for d_atom in atom:
                d_atom.coord = coord_operation(d_atom.coord, key, op_n, distance =distance)
        atom.coord = coord_operation(atom.coord, key, op_n, distance=distance)
    return entity


def convertFromOrthToFrac(orth_coords:list[float], parameters:dict) -> list:
    """
    Convert orthogonal coordinate to fractional coordinate.
    :param orth_coords: Orthogonal coordinate
    :param parameters: Operation parameters
    :return: Fractional coordinate
    """
    assert len(orth_coords) == 3
    x, y, z = orth_coords

    nx = (x * parameters["vvy"]) + (y * parameters["vvz"]) + (z * parameters["uuz"])
    ny = (y * parameters["uuy"]) + (z * parameters["vv"])
    nz = z * parameters["uu"]

    return [nx, ny, nz]


def convertFromFracToOrth(frac_coords:list[float], parameters:dict) -> list:
    """
    Convert fractional coordinate to orthogonal coordinate.
    :param frac_coords: Fractional Coordinate
    :param parameters: Operation parameters
    :return: Orthogonal coordinate
    """
    assert len(frac_coords) == 3
    t1, t2, t3 = frac_coords

    tz = t3 / parameters["uu"]
    ty = (t2 - tz * parameters["vv"]) / parameters["uuy"]
    tx = (t1 - ty * parameters["vvz"] - tz * parameters["uuz"]) / parameters["vvy"]

    return tx, ty, tz


def entity_to_frac(entity:bp.Entity.Entity, parameters:dict) -> bp.Entity.Entity:
    """
    Convert orthogonal Entity coordinates to fractional Entity coordinates in-place.
    :param entity: Original orthogonal Entity
    :param parameters: Operation parameters
    :return: Original Entity with fractional coordinates.
    """
    for atom in entity.get_atoms():
        if atom.is_disordered() > 0:
            for d_atom in atom:
                d_atom.coord = convertFromOrthToFrac(d_atom.coord, parameters)
        atom.coord = convertFromOrthToFrac(atom.coord, parameters)
    return entity


def entity_to_orth(entity:bp.Entity.Entity, parameters:dict) -> bp.Entity.Entity:
    """
    Convert fractional Entity coordinates to orthogonal Entity coordinates in-place.
    :param entity: Original fractional Entity
    :param parameters: Operation parameters
    :return: Original Entity with orthogonal coordinates.
    """
    for atom in entity.get_atoms():
        if atom.is_disordered() > 0:
            for d_atom in atom:
                d_atom.coord = convertFromFracToOrth(d_atom.coord, parameters)
        atom.coord = convertFromFracToOrth(atom.coord, parameters)
    return entity


def generate_displaced_copy(original:bp.Entity.Entity, distance:list[float]|float = 99.5, key:int = None, op_n:int = None) -> bp.Entity.Entity:
    """
    Create a copy of an Entity displaced a defined distance in-place.
    :param original: Original Entity
    :param distance: Distance to displace entity, float or [X, Y, Z]
    :param key: Space group key
    :param op_n: Operation number
    :return: Original Entity displaced
    """
    if original is None:
        return None
    displaced = original.copy()

    if distance is None:
        distance = [0, 0, 0]

    if type(distance) is float or type(distance) is int:
        distance = [distance, distance, distance]

    if key is None and op_n is None:
        for atom in displaced.get_atoms():
            if atom.is_disordered() > 0:
                for d_atom in atom:
                    d_atom.coord = [x+d for x, d in zip(d_atom.coord, distance)]
            else:
                atom.coord = [x+d for x, d in zip(atom.coord, distance)]
    else:
        coord_operation_entity(displaced, key=key, op_n=op_n, distance =distance)
        for atom in displaced.get_atoms():
            if atom.is_disordered() > 0:
                pass

    return displaced




def get_fractional_distance(coord1, coord2, params, root=False) -> float:
    """
    Calculate the square distance in orthogonal space between two coordinates in fractional space.
    Optionally returns the root square distance.
    :param coord1: First coordinate.
    :param coord2: Second coordinate.
    :param params: Crystal parameters.
    :param root: Whether to return the root square distance.
    :return: Orthogonal square distance (or root square distance).
    """
    deltaX, deltaY, deltaZ = vector(coord1, coord2)

    a = params["A"]
    b = params["B"]
    c = params["C"]
    c_a = params["c_a"]
    c_b = params["c_b"]
    c_g = params["c_g"]

    d2 = (a**2)*(deltaX**2) + (b**2)*(deltaY**2) + (c**2)*(deltaZ**2) +2*b*c*c_a*deltaY*deltaZ +2*a*c*c_b*deltaX*deltaZ +2*a*b*c_g*deltaX*deltaY
    if root:
        return np.sqrt(*d2)
    else:
        return d2



