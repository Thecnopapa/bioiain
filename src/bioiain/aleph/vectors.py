import os, json
import numpy as np

from ..utilities.logging import log
from ..utilities.exceptions import *
from ..utilities.maths import *
from ..visualisation.plots import plot_heatmap








class CVector(object):
    def __init__(self, res1, res2, res3):
        self.res1 = res1
        self.res2 = res2
        self.res3 = res3
        self.residues = (res1, res2, res3)

        self.chain = res2.chain
        self.resnum = res2.resnum

        self.id = (self.resnum, self.chain)


        self.d = None
        self.start = None
        self.end = None

        self.calculate()


    def __repr__(self):
        return f"<bi.{self.__class__.__name__} c.{self.chain} ({self.res1.resnum}>{self.res2.resnum}>{self.res3.resnum})>"


    def _centroid(self, *coords):
        coords = np.array(coords)
        coords = coords.sum(axis=0)
        coords = coords/len(coords)
        return coords



    def calculate(self):

        self.start = find_com([r.ca for r in self.residues])
        self.end = find_com([r.o for r in self.residues])

        self.v = vector(self.start, self.end)


        self.d = length(self.v)


        return self

    def __mod__(self, target:CVector):
        return CVPair(self, target)




class CVPair(object):
    def __init__(self, cvec1, cvec2):
        self.v1 = cvec1
        self.v2 = cvec2

        self.chain1 = self.v1.chain
        self.chain2 = self.v2.chain

        self.resnum1 = self.v1.resnum
        self.resnum2 = self.v2.resnum


        self.v = None
        self.d = None
        self.a = None
        self.t1 = None
        self.t2 = None

        self.calculate()

    def __repr__(self):
        return f"<bi.{self.__class__.__name__} {self.v1.id} - {self.v2.id}>"


    def calculate(self):


        self.v = vector(self.v1.start, self.v2.start)
        self.d = length(self.v)
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
        for n, v1 in enumerate(self.vectors[:-1]):
            self.matrix.append([])
            for v2 in self.vectors[n+1:]:
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


    def save_fig(self, attribute="d", save_folder="./cvmaps"):
        os.makedirs(save_folder, exist_ok=True)
        filename = os.path.join(save_folder, f"{attribute}.png")
        plot_heatmap(self.square(attribute), filename=filename)
        return filename

    def show(self, attribute="d"):
        plot_heatmap(self.square(attribute), show=True)

    def calculate_neighbours(self):
        for n1, cv in enumerate(self.vectors):
            print("NVECTORS:", len(self.vectors))

            target = (99999., None)
            for n2, vp in enumerate(self.matrix[n1]):


                if vp is None:
                    vp = self.matrix[n2][n1]

                #print(vp.chain1, vp.chain2)

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
                    target = (vp.d, t)
                    continue

               
            cv.closest = target
            print("closest to", cv, "is", cv.closest)
            if target[1] is None:
                raise Exception("AAAAA")






