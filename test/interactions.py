### src.bioiain import #################################################################################################
import sys
sys.path.append('..')
from src.bioiain import *
from src.bioiain.utilities import *
########################################################################################################################




class PPI(object):
    def __init__(self, **kwargs):
        self.atoms1:list|None = None
        self.atoms2:list|None = None

        self.fragments1:list|None = None
        self.fragments2:list|None = None

    def __repr__(self):
        return f"<bi.PPI at {hex(id(self))}>"

    @classmethod
    def from_atoms(cls, atoms1:list, atoms2:list, **kwargs):
        self = cls(**kwargs)
        self.atoms1 = atoms1
        self.atoms2 = atoms2
        return self

    @classmethod
    def from_fragments(cls, fragments1:list, fragments2:list, **kwargs):
        self = cls(**kwargs)
        self.fragments1 = fragments1
        self.fragments2 = fragments2

        self.atoms1 = []
        self.atoms2 = []

        for frag in self.fragments1:
            self.atoms1.extend(frag.atoms())
        for frag in self.fragments2:
            self.atoms2.extend(frag.atoms())  
        return self



def get_all_PPIs(entity):
    fragment_interactions = []
    kdtree = entity.ca_kdtree()
    print(kdtree)
    for res in entity.residues():
        ca = res.ca
        nns = kdtree.radius(ca)
        f1 = ca.fragment()
        for nn in nns[0]:
            #print(nn)
            natom = kdtree.atom_of(nn)
            #print(natom.fragment())
            f2 = natom.fragment()
            fragment_interactions.append([f1,f2])
    print(fragment_interactions)




if __name__ == "__main__":

    entity = base.entity.Entity.from_file("./1M2Z.cif")
    entity.fragment(in_place=True)
    print(entity)

    get_all_PPIs(entity)




