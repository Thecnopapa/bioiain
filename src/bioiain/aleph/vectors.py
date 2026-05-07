import os, json
import numpy as np

from .. import SUBDIR_NAME
from ..utilities.logging import log
from ..utilities.exceptions import *
from ..utilities.maths import *
from ..visualisation.plots import plot_heatmap
from ..base.atom import PseudoAtom
from ..tools.SASA import KDT







class CVector(object):
    def __init__(self, res1, res2, res3, params=None, symops=None, entity_centre=None, vc_mode=None):
        self.res1 = res1
        self.res2 = res2
        self.res3 = res3
        self.residues = (res1, res2, res3)

        self.chain = res2.chain
        self.resnum = res2.resnum
        self.resname = res2.resname
        self.resseq = res2.resseq
        self.fragment = res2.fragment

        self.id = (self.resnum, self.chain)
        self.across_chains = False
        self.is_gap = False
        self.trash = False

        self.vc_mode = vc_mode
        self.vc = None
        self.d = None
        self.v = None
        self.start = None
        self.end = None

        self.closest = None
        self.closest_vp = None
        self.closest_opn = None
        self.closest_pos = None

        self.closest_lig = None
        self.dist_to_lig = None

        self.params = params
        self.symops = symops
        self.entity_centre = entity_centre

        self.calculate()

    def full_id(self):
        return f"{self.chain}_{self.fragment}_{self.resseq}_{self.resname}_{self.resnum}"

    def _mmcif_dict(self):
        d =  dict(
            id=self.full_id(),
            chain=self.chain,
            fragment=self.fragment,
            resseq=self.resseq,
            resname=self.resname,
            resnum=self.resnum,
            d=f"{self.d:8.3f}",
            start_x=f"{self.start.coord[0]:8.3f}", start_y=f"{self.start.coord[1]:8.3f}", start_z=f"{self.start.coord[2]:8.3f}",
            end_x=f"{self.end.coord[0]:8.3f}", end_y=f"{self.end.coord[1]:8.3f}", end_z=f"{self.end.coord[2]:8.3f}",
            vc_x=f"{self.vc.coord[0]:8.3f}", vc_y=f"{self.vc.coord[1]:8.3f}", vc_z=f"{self.vc.coord[2]:8.3f}",
            closest_resseq=self.closest.resseq if self.closest is not None else None,
            closest_resnum=self.closest.resnum if self.closest is not None else None,
            closest_chain=self.closest.chain if self.closest is not None else None,
            closest_fragment=self.closest.fragment if self.closest is not None else None,
            closest_pos="_".join([str(p) for p in self.closest_pos]) if self.closest_pos is not None else None,
            closest_opn=self.closest_opn,
            closest_lig_name=self.closest_lig.name if self.closest_lig is not None else None,
            closest_lig_chain=self.closest_lig.chain if self.closest_lig is not None else None,
            dist_to_lig=f"{self.dist_to_lig:8.3f}" if self.dist_to_lig is not None else None,
        )
        for k, v in d.items():
            if v is None:
                d[k] = "."
            d[k] = str(d[k])
        return d


    def __repr__(self):
        return f"<bi.{self.__class__.__name__} c.{self.chain} ({self.res1.resnum}>{self.res2.resnum}>{self.res3.resnum})>"


    def _centroid(self, *coords):
        coords = np.array(coords)
        coords = coords.sum(axis=0)
        coords = coords/len(coords)
        return coords



    def calculate(self):

        self.start =  PseudoAtom(find_com([r.ca for r in self.residues]))
        self.end = PseudoAtom(find_com([r.o for r in self.residues]))

        self.v = vector(self.start.coord, self.end.coord)
        self.d = length(self.v)
        if self.vc_mode is None:
            self.vc = self.start.copy()
        elif self.vc_mode == "ca_projection":
            ca = self.res2.ca.coord
            start_ca = vector(self.start.coord, ca)
            self.vc =  PseudoAtom((start_ca[0] + ca[0], start_ca[1] + ca[1], start_ca[2] + ca[2]))
        else:
            log("error", "VC Mode Not Implemented:", self.vc_mode)
            raise NotImplementedError("VC Mode Not Implemented:", self.vc_mode)

        if len({self.res1.chain, self.res2.chain, self.res3.chain}) > 1:
            self.across_chains = True
            self.trash = True

        if abs(self.res1.resnum - self.res2.resnum) != 1 or abs(self.res2.resnum - self.res3.resnum) != 1:
            self.is_gap = True
            self.trash = True


        return self

    def __mod__(self, target):
        return CVPair(self, target)




class CVPair(object):
    def __init__(self, cvec1, cvec2):
        self.v1 = cvec1
        self.v2 = cvec2
        self.opn_of_v2 = None
        self.pos_of_v2 = None

        self.chain1 = self.v1.chain
        self.chain2 = self.v2.chain

        self.resnum1 = self.v1.resnum
        self.resnum2 = self.v2.resnum

        self.fragment1 = self.v1.fragment
        self.fragment2 = self.v2.fragment

        self.dlig = None

        self.v = None
        self.d = None
        self.a = None
        self.t1 = None
        self.t2 = None

        self.calculate()

    def __repr__(self):
        if self.fragment1 is not None and self.fragment2 is not None:
            return f"<bi.{self.__class__.__name__} {self.v1.resnum}(F{self.v1.fragment}) - {self.v2.resnum}(F{self.v2.fragment})>"
        else:
            return f"<bi.{self.__class__.__name__} {self.v1.resnum}({self.v1.chain}) - {self.v2.resnum}({self.v2.chain})>"


    def calculate(self):


        
        if self.v2.params is not None and self.v2.symops is not None:
            self.v1.vc.to_frac(self.v2.params)
            self.v2.vc.to_frac(self.v2.params)
            coord, d, opn, pos = self.v2.vc.closest(self.v1.vc.coord, symops=self.v2.symops, params=self.v2.params, centre=self.v2.entity_centre)

            self.d = d
            self.opn_of_v2 = opn
            self.pos_of_v2 = pos
            self.v1.vc.to_orth(self.v2.params)
            self.v2.vc.to_orth(self.v2.params)
            self.v = vector(self.v1.vc.coord, coord)
        else:
            self.v = vector(self.v1.vc.coord, self.v2.vc.coord)
            self.d = flength(self.v)


        #print(self.v, self.d)
        self.a = angle_between_vectors(self.v1.v, self.v2.v)
        self.t1 = angle_between_vectors(self.v1.v, self.v)
        self.t2 = angle_between_vectors(self.v2.v, self.v)
        self.da = dihedral_angle(self.v1.end, self.v1.start, self.v2.start, self.v2.end)
        return self

    def map_lig(self):
        try:
            self.dlig = min([d for d in (self.v1.dist_to_lig, self.v2.dist_to_lig) if d is not None])
        except:
            raise

    def _mmcif_dict(self):
        data = {
            "cv1_id": self.v1.full_id(),
            "cv2_id": self.v2.full_id(),
            "resnum_1": self.v1.resnum,
            "chain1": self.v1.chain,
            "fragment1": self.v1.fragment,
            "resnum_2": self.v2.resnum,
            "chain2": self.v2.chain,
            "fragment2": self.v2.fragment,
            "opn_of_v2": self.opn_of_v2,
            "pos_of_v2": self.v2.pos,
            "d": f"{self.d:8.3f}",
            "a": f"{self.a:8.3f}",
            "t1": f"{self.t1:8.3f}",
            "t2": f"{self.t2:8.3f}",
            "da": f"{self.da:8.3f}",
            "dlig": f"{self.dlig:8.3f}",
            "v_x": f"{self.v[0]:8.3f}",
            "v_y": f"{self.v[1]:8.3f}",
            "v_z": f"{self.v[2]:8.3f}"
        }
        return data




class CVMatrix(object):
    def __init__(self, cvector_list, vc_mode=None, complete=False, max_distance=30, **kwargs):
        self.vectors = cvector_list
        self.matrix = None
        self.length = len(self.vectors)
        self.vc_mode = vc_mode
        self.max_distance = max_distance

        self.calculate_neighbours(max_distance=max_distance, **kwargs)


    def calculate_full(self):
        self.matrix = []
        log(2, f"Calculating matrix ({self.vc_mode})...")

        for n, v1 in enumerate(self.vectors):
            self.matrix.append([])
            print(f"{n:5d}/{len(self.vectors):5d}", end="\r")

            for v2 in self.vectors:
                if v1 == v2:
                    pair = None
                else:
                    pair = v1 % v2
                self.matrix[n].append(pair)

        for n, r in enumerate(self.matrix):
            extension = self.length - len(r) -1
            zeros = [None] * extension
            self.matrix[n] = zeros + self.matrix[n]

        log(3, "n vectors", len(self.vectors))
        log(3 , "m length", len(self.matrix))
        return self

    def map(self, attribute):
        log(2, f"Mapping CVMatrix (attr={attribute})")
        return (lambda m: np.array([np.array([getattr(x, attribute) if x is not None and getattr(x, attribute) is not None else 0 for x in l]) for l in m ]))(self.matrix)

    def square(self, attribute):

        t = self.map(attribute)
        sq = square_matrix(t)
        return sq

    def triangle(self, attribute="d"):
        t = self.map(attribute)
        return t


    def save_fig(self, attribute="d", save_folder=None, prefix="cvmatrix"):
        if save_folder is None:
            save_folder = os.path.join(SUBDIR_NAME, "cvmaps")
        os.makedirs(save_folder, exist_ok=True)
        if self.vc_mode is not None:
            prefix = f"{prefix}_{self.vc_mode}"
        filename = os.path.join(save_folder, f"{prefix}_{attribute}.png")
        plot_heatmap(self.square(attribute), filename=filename)
        return filename

    def show(self, attribute="d"):
        plot_heatmap(self.square(attribute), show=True)

    def map_ligands(self, entity):
        log(2, "Mapping ligands...")
        ligands = entity.ligands()
        if len(ligands) == 0:
            log("warning", f"No ligands found in {entity.name()}")
            for cv in self.vectors:
                cv.dist_to_lig = 99
                cv.closest_lig = None
            return self

        for cv in self.vectors:
            distances = [length(vector(cv.start, l.com())) for l in ligands]
            #print(distances)
            cv.dist_to_lig = min(distances)
            cv.closest_lig = ligands[distances.index(cv.dist_to_lig)]
            #log(1, f"Closest lig ({cv}): {cv.closest_lig} at {cv.dist_to_lig:4.2f}A")

        for row in self.matrix:
            for vp in row:
                if vp is None:
                    continue
                vp.map_lig()


        return self



    def calculate_neighbours(self, use_fragments=True, max_distance=30, n_neighbours=1):
        log(2, "Calculating neighbours for:", self)


        data = {k: {"vc":cv.vc, "cv":cv} for k, cv in enumerate(self.vectors)}
        for k, v in data.items():
            v["vc"]._nnk = k


        tree = KDT([v["vc"] for v in data.values()])

        for n1, v in data.items():
            cv = v["cv"]
            vc = v["vc"]


            radius = 5
            neighs_found=False
            while radius <= max_distance:
                neighbors, distances = tree.of(vc, radius=radius, distances=True)
                print(neighbors)
                print(distances)
                if len(neighbors) == 0:
                    log("error", f"No neighbours found radius={radius}, even itself")
                    raise Exception()
                if len(neighbors) < n_neighbours+1:
                    radius += 5
                    continue

                neighs_found = True
                targets = [(tree.atom_of(neigh), d) for neigh, d in zip(neighbors, distances)]
                targets = sorted(targets, key=lambda x: x[0])
                print("sorted targets:")
                print(targets)

                for t, d in targets[1:n_neighbours+1]:
                    k = t._nnk
                    cv.closest = data[k]["cv"]
                    cv.closest_vp = cv % cv.closest
                    cv.closest_opn = cv.closest_vp.opn_of_v2
                    cv.closest_pos = cv.closest_vp.pos_of_v2
                    if n_neighbours != 1:
                        log("Warning", "More than one neighbour requested but not implemented")
                    break # Only programmed for n_neighbours == 1, must be adapted if more needed
                break

            if not neighs_found:
                log("ERROR", f"No neighbours found for cv: {cv}\nMight be due to small structure")
                raise NoNeighboursFound()
            return self








            for n2, vp in enumerate(self.matrix[n1]):

                if vp is None:
                    vp = self.matrix[n2][n1]

                if n2 == n1 or (vp.chain1 == vp.chain2 and vp.opn_of_v2 is not None and vp.pos_of_v2 is not None and abs(vp.resnum1-vp.resnum2) <=2) :
                    continue
                print(f"{n1+1} / {n2+1} / {len(self.vectors)}", end="\r")


                    #if vp is None:
                    #    continue
                #print(vp, vp.d)
                #print(vp.chain1, vp.chain2)

                if use_fragments:
                    #print(vp.fragment1, vp.fragment2)
                    if vp.fragment1 == vp.fragment2:
                        continue
                else:
                    if vp.chain1 == vp.chain2:
                        continue


                #print(str(cv), str(vp.v1), str(vp.v2))

                if str(cv) == str(vp.v1):
                    t = vp.v2
                elif str(cv) == str(vp.v2):
                    t = vp.v1
                else:
                    continue

                #print(cv, t)
                if vp.d < target[0]:
                    target = (vp.d, t, vp, vp.opn_of_v2, vp.pos_of_v2)
                    continue

            #print(target)
            cv.closest = target[1]
            cv.closest_vp = target[2]
            cv.closest_opn = target[3]
            cv.closest_pos = target[4]
            #print("closest to", cv, "is", cv.closest)
            if target[1] is None:
                log("ERROR", f"No neighbours found for cv: {cv}\nMight be due to small structure")
                raise NoNeighboursFound()

        return self






