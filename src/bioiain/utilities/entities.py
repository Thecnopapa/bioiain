import Bio.PDB as bp



def print_all_coords(entity:bp.Entity.Entity, head:int = 5) -> None:
    """
    Print all coordinates of an entity.
    :param entity: Entity to print
    :param head: N coordinates to print (default 5), -1 to print all
    """
    n = 0
    for atom in entity.get_atoms():
        print(atom.coord)
        n += 1
        if n >= head:
            break
