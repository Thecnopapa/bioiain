import os, json
import numpy as np

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
        self.fragment = res2.fragment

        self.id = (self.resnum, self.chain)


        self.d = None
        self.start = None
        self.end = None

        self.closest = None
        self.closest_vp = None

        self.params = params
        self.symops = symops
        self.entity_centre = entity_centre

        self.calculate()


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




class CVMatrix(object):
    def __init__(self, cvector_list):
        self.vectors = cvector_list
        self.matrix = None
        self.length = len(self.vectors)

        self.calculate()

    def calculate(self):
        self.matrix = []
        print("CALCULATING matrix...")
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


        return self

    def map(self, attribute):
        return (lambda m: np.array([np.array([getattr(x, attribute) if x is not None else 0 for x in l]) for l in m ]))(self.matrix)

    def square(self, attribute):

        t = self.map(attribute)
        sq = square_matrix(t)
        return sq

    def triangle(self, attribute="d"):
        t = self.map(attribute)
        return t


    def save_fig(self, attribute="d", save_folder="./bioiain/cvmaps"):
        os.makedirs(save_folder, exist_ok=True)
        filename = os.path.join(save_folder, f"{attribute}.png")
        plot_heatmap(self.square(attribute), filename=filename)
        return filename

    def show(self, attribute="d"):
        plot_heatmap(self.square(attribute), show=True)

    def calculate_neighbours(self, use_fragments=True):
        print("MAPPING neighbours")
        print("n vectors", len(self.vectors))
        print("m shape", len(self.matrix))

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






