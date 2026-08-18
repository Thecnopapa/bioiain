### src.bioiain import #################################################################################################
import sys
sys.path.append('..')
from src.bioiain import *
from src.bioiain.utilities import *
########################################################################################################################

from src.bioiain.base import Entity
import networkx as nx
from src.bioiain.visualisation.plots import fig2D, show, close

class PPI(object):
    def __init__(self, **kwargs):
        self.atoms1:list|None = None
        self.atoms2:list|None = None

        self.fragments1:list|None = None
        self.fragments2:list|None = None

        self.subentity1:Entity|None = None
        self.subentity2:Entity|None = None

        self.op: int|None = kwargs.get("op", None)
        self.name: str = kwargs.get("name", "")

    def __repr__(self):
        if self.atoms1 is None or self.atoms2 is None:
            return f"<bi.PPI at {hex(id(self))}>"
        else:
            return f"<bi.PPI {self.name} {len(self.atoms1)}x{len(self.atoms2)} atoms (op:{self.op})>"

    def _init_subentities(self):
        self.subentity1 = Entity.from_atoms(self.atoms1, code=self.name.split("_")[0], name=f"{self.name}_p1", share=False)
        self.subentity2 = Entity.from_atoms(self.atoms2, code=self.name.split("_")[0], name=f"{self.name}_p2", share=False)

        self.subentity1.paths["sub_folder"] = "fragmented/PPIs"
        self.subentity2.paths["sub_folder"] = "fragmented/PPIs"

        self.subentity1.extension = "ppi"
        self.subentity2.extension = "ppi"

        self.subentity1.export()
        self.subentity2.export()
        return self

    @classmethod
    def from_atoms(cls, atoms1:list, atoms2:list, **kwargs):
        self = cls(**kwargs)
        self.atoms1 = atoms1
        self.atoms2 = atoms2
        self._init_subentities()
        
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
        self._init_subentities()
        return self



def get_all_PPIs(entity, radius=10, min_contacts=1, intra_asu=True):
    fragment_interactions = []
    interaction_dict = {}
    ppis = []
    kdtree = entity.ca_kdtree()
    #print(kdtree)
    for res in entity.residues():
        ca = res.ca
        nns = kdtree.radius(ca, radius=radius)
        f1 = ca.fragment()
        for nn in nns[0]:
            #print(nn)
            natom = kdtree.atom_of(nn)
            #print(natom.fragment())
            f2 = natom.fragment()
            op2 = kdtree.op_of(nn)
            pos2 = kdtree.pos_of(nn)
            is_intra=False
            if op2 == 1:
                is_intra=True
                if intra_asu:
                    if ca.chain == natom.chain:
                        continue
                    if f2 <= f1:
                        continue
                else:
                    continue
            key = (f1,f2, pos2)
            fragment_interactions.append([f1,f2, op2, pos2])
            if op2 not in interaction_dict:
                interaction_dict[op2] = {}
            if key not in interaction_dict[op2]:
                interaction_dict[op2][key] = {"n":0, "c1":ca.chain, "c2": natom.chain, "intra":is_intra, "f1":f1, "f2":f2, "op": op2, "pos": pos2}
            interaction_dict[op2][key]["n"] += 1

    for o in interaction_dict.keys():
        interaction_dict[o] = {k:v for k,v in sorted(sorted(interaction_dict[o].items(), key=lambda x: x[0][1]), key=lambda x: x[0][0]) if v["n"] >= min_contacts}
    interaction_dict = {k:v for k,v in sorted(interaction_dict.items(), key = lambda x: x[0]) if len(v) > 0}
    #print(fragment_interactions)
    #print(interaction_dict)
    
    for o, vv in interaction_dict.items():
        graph = nx.Graph()
        
        #graph.add_node(f"op{o}", color="red")
        for k, v in vv.items():
            k1 = f"{v['c1']}{v["f1"]}"
            # if v["intra"]:
            #     k2 = f"s{v['c2']}{v["f2"]}{v["pos"]}"
            # else:
            k2 = f"s{v['c2']}{v["f2"]}{v["pos"]}"
            graph.add_node(k1, color="blue")
            graph.add_node(k2, color="green")
            graph.add_edge(k1, k2, weight=v["n"])
            #graph.add_edge(k[0], f"op{o}", weight=1)
        #print(graph)
        for n, cc in enumerate(nx.connected_components(graph)):
            fig, ax = fig2D(figsize=(10,10))
            subgraph = graph.subgraph(cc)
            #print(subgraph)
            #print(subgraph.nodes)
        
            nx.draw(subgraph, with_labels=True, font_weight='bold', font_size=8, ax=ax, node_color=dict(subgraph.nodes.data("color")).values())
            ax.set_aspect('equal')
            ax.set_title(f"ppi:{o}-{n}")
            os.makedirs(os.path.join(TEMP_FOLDER,"graphs"), exist_ok=True)
            fig_path = os.path.join(TEMP_FOLDER,"graphs", f"graph{o}-{n}.png")
            fig.savefig(fig_path)
            print("open", fig_path)

            atoms1 = []
            atoms2 = []
            print(cc)
            for atom, op, pos in zip(kdtree.atoms, kdtree.operations, kdtree.positions):
                if (op == 1) and f"{atom.chain}{atom.fragment()}" in cc:
                    atoms1.extend(atom._residue.atoms)
                elif (op == o) and f"s{atom.chain}{atom.fragment()}{pos}" in cc:
                    atoms2.extend(a.copy().symop(symop=entity.symops(op), params=entity.params(), position=pos) for a in atom._residue.atoms)
            #print(atoms1)
            #print(atoms2)

            ppi = PPI.from_atoms(atoms1, atoms2, op=o, name=f"{entity.name()}_{o}-{n}")
            ppis.append(ppi)
            print(ppi)
    return ppis


def plot_ppis(entity, ppis):
    from src.bioiain.visualisation.pymol import PymolScript
    script = PymolScript(name=f"{entity.name()}_ppis", use_temp=True)
    print(script)
    script.load(entity.path())
    for ppi in ppis:
        print(ppi)
        ppi.subentity1.path()
        script.load(ppi.subentity1.path())
        script.load(ppi.subentity2.path())
        script.group(ppi.name, ppi.name)
    script.orient()
    script.execute()





if __name__ == "__main__":

    entity = base.entity.Entity.from_file("./1M2Z.cif")
    entity.fragment(in_place=True)
    print(entity)

    ppis = get_all_PPIs(entity, min_contacts=2)

    plot_ppis(entity, ppis)




