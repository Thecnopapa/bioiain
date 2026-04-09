import os, json
import numpy as np

from .. import SUBDIR_NAME
from ..utilities.logging import log
from ..utilities.exceptions import *
from ..utilities.maths import *
from ..visualisation.plots import plot_heatmap
from ..base.atom import PseudoAtom







class CVector(object):
    def __init__(self, res1, res2, res3, params=None, symops=None, entity_centre=None):
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


    def _mmcif_dict(self):
        d =  dict(
            chain=self.chain,
            fragment=self.fragment,
            resseq=self.resseq,
            resname=self.resname,
            resnum=self.resnum,
            d=f"{self.d:8.3f}",
            start_x=f"{self.start.coord[0]:8.3f}", start_y=f"{self.start.coord[1]:8.3f}", start_z=f"{self.start.coord[2]:8.3f}",
            end_x=f"{self.end.coord[0]:8.3f}", end_y=f"{self.end.coord[1]:8.3f}", end_z=f"{self.end.coord[2]:8.3f}",
            closest_resseq=self.closest.resseq if self.closest is not None else None,
            closest_resnum=self.closest.resnum if self.closest is not None else None,
            closest_chain=self.closest.chain if self.closest is not None else None,
            closest_fragment=self.closest.fragment if self.closest is not None else None,
            closest_pos="_".join([str(p) for p in self.closest_pos]) if self.closest_pos is not None else None,
            closest_opn=self.closest_opn,
            closest_lig_name=self.closest_lig.name if self.closest_lig is not None else None,
            closest_lig_chain=self.closest_lig.chain if self.closest_lig is not None else None,
            dist_to_lig=f"{self.dist_to_lig:8.3f}",
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
            self.v1.start.to_frac(self.v2.params)
            self.v2.start.to_frac(self.v2.params)
            coord, d, opn, pos = self.v2.start.closest(self.v1.start.coord, symops=self.v2.symops, params=self.v2.params, centre=self.v2.entity_centre)

            self.d = d
            self.opn_of_v2 = opn
            self.pos_of_v2 = pos
            self.v1.start.to_orth(self.v2.params)
            self.v2.start.to_orth(self.v2.params)
            self.v = vector(self.v1.start.coord, coord)
        else:
            self.v = vector(self.v1.start.coord, self.v2.start.coord)
            self.d = flength(self.v)


        #print(self.v, self.d)
        self.a = angle_between_vectors(self.v1.v, self.v2.v)
        self.t1 = angle_between_vectors(self.v1.v, self.v)
        self.t2 = angle_between_vectors(self.v2.v, self.v)
        return self

    def map_lig(self):
        try:
            self.dlig = min([d for d in (self.v1.dist_to_lig, self.v2.dist_to_lig) if d is not None])
        except:
            raise

class CVMatrix(object):
    def __init__(self, cvector_list):
        self.vectors = cvector_list
        self.matrix = None
        self.length = len(self.vectors)

        self.calculate()

    def calculate(self):
        self.matrix = []
        log(2, "Calculating matrix...")

        for n, v1 in enumerate(self.vectors):
            self.matrix.append([])
            print(f"{n} / {len(self.vectors)}", end="\r")

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


    def save_fig(self, attribute="d", save_folder=None):
        if save_folder is None:
            save_folder = os.path.join(SUBDIR_NAME, "cvmaps")
        os.makedirs(save_folder, exist_ok=True)
        filename = os.path.join(save_folder, f"{attribute}.png")
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



    def calculate_neighbours(self, use_fragments=True):
        log(2, "Calculating neighbours for:", self)


        for n1, cv in enumerate(self.vectors):

            target = (99999., None, None, None, None)
            for n2, vp in enumerate(self.matrix[n1]):

                if vp is None:
                    vp = self.matrix[n2][n1]

                if n2 == n1 or (vp.chain1 == vp.chain2 and vp.opn_of_v2 is not None and vp.pos_of_v2 is not None and abs(vp.resnum1-vp.resnum2) <=2 ) :
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
                raise Exception("AAAAA")

        return self






