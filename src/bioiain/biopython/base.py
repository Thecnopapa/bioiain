


import Bio.PDB as bp



class BiopythonOverlayClass:
    @classmethod
    def cast(cls, entity:bp.Entity.Entity):
        entity.__class__ = cls
        return entity

