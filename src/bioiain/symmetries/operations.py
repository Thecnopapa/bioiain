import numpy as np
from space_groups import dictio_space_groups




def convertFromOrthToFrac(orth_coords, parameters):
    x, y, z = orth_coords

    nx = (x * parameters["vvy"]) + (y * parameters["vvz"]) + (z * parameters["uuz"])
    ny = (y * parameters["uuy"]) + (z * parameters["vv"])
    nz = z * parameters["uu"]

    return nx, ny, nz


def entity_to_frac(entity, params):
    for atom in entity.get_atoms():
        if atom.is_disordered() > 0:
            for d_atom in atom:
                d_atom.coord = convertFromOrthToFrac(d_atom.coord, params)
        #else:
        atom.coord = convertFromOrthToFrac(atom.coord, params)
    return entity


def convertFromFracToOrth(frac_coords, parameters):
    t1, t2, t3 = frac_coords

    tz = t3 / parameters["uu"]
    ty = (t2 - tz * parameters["vv"]) / parameters["uuy"]
    tx = (t1 - ty * parameters["vvz"] - tz * parameters["uuz"]) / parameters["vvy"]

    return tx, ty, tz

def entity_to_orth(entity, params):
    for atom in entity.get_atoms():
        if atom.is_disordered() > 0:
            for d_atom in atom:
                d_atom.coord = convertFromFracToOrth(d_atom.coord, params)
        #else:
        atom.coord = convertFromFracToOrth(atom.coord, params)
    return entity


def generate_displaced_copy(original, distance = 99.5, key = None, op_n = None):
    #print6("Generating displaced copy")
    if original is None:
        return None
    displaced = original.copy()

    try:
        _ = [x for x in distance]
    except:
        distance = [distance] * 3
    #print6(distance)

    if key is None and op_n in None:
        for atom in displaced.get_atoms():
            if atom.is_disordered() > 0:
                for d_atom in atom:
                    d_atom.coord = [x+d for x, d in zip(d_atom.coord, distance)]
            else:
                atom.coord = [x+d for x, d in zip(atom.coord, distance)]
    else:
        #print(key,op_n)
        coord_operation_entity(displaced, key=key, op_n=op_n, distance =distance)
        for atom in displaced.get_atoms():
            if atom.is_disordered() > 0:
                #print(atom.coord)
                #[print(d_atom.coord) for d_atom in atom]
                pass

    return displaced


def coord_operation_entity(entity, key, op_n, distance = (0,0,0)):
    for atom in entity.get_atoms():
        if atom.is_disordered() > 0:
            #print(atom.is_disordered())
            for d_atom in atom:
                #print(d_atom.get_full_id())
                d_atom.coord = coord_operation(d_atom.coord, key, op_n, distance =distance)
        #else:
        atom.coord = coord_operation(atom.coord, key, op_n, distance=distance)
    return entity

def coord_operation(coord, key, op_n, distance = (0,0,0)):

    rotation = dictio_space_groups[key]
    #print(rotation)
    operation = rotation["symops"][op_n]
    rot = operation["rot"]
    tra = operation["tra"]

    x, y, z = coord

    nx = (rot[0][0] * x) + (rot[0][1] * y) + (rot[0][2] * z) + tra[0]+distance[0]
    ny = (rot[1][0] * x) + (rot[1][1] * y) + (rot[1][2] * z) + tra[1]+distance[1]
    nz = (rot[2][0] * x) + (rot[2][1] * y) + (rot[2][2] * z) + tra[2]+distance[2]

    return nx, ny, nz

def coord_add(coord, deltas, subtract = False):
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

def get_operation(key,op_n):
    return dictio_space_groups[key]["symops"][op_n]