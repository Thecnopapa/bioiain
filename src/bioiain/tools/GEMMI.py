import gemmi







def gemmi_cell(cell ,space_group):
    if type(cell) in (list, tuple):
        cell = gemmi.UnitCell(*cell)
    if type(space_group) is str:
        space_group = gemmi.SpaceGroup(space_group) # 'I 2 2 2'
    gv = gemmi.GruberVector(cell, space_group)
    return gv



def niggli_cell(gv):
    gv.niggli_reduce()
    return gv.get_cell()

def buerger_cell(gv):
    gv.buerger_reduce()
    return gv.get_cell()
