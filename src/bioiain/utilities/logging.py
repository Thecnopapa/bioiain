def print_all_coords(entity, head = 5):
    n = 0
    for atom in entity.get_atoms():
        print(atom.coord)
        n += 1
        if n >= head:
            break


