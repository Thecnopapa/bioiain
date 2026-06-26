import os, sys, json
from ..base import *
from ..utilities.exceptions import *
from ..utilities.logging import log

import numpy as np


class _ALEPH(object):

    def __init__(self, params={}):
        self.params = dict(
            min_ah=4, min_bs=3,
            strictness_ah=0.45, strictness_bs=0.20,
            weight = "distance_avg",
        )
        self.params |= params

        self.dist_mean_bs = 5.1
        self.dist_mean_ah = 0.0 # Unused in ALEPH
        self.dist_num_ah = 0 # Unused in ALEPH
        self.angle_mean_bs = 54
        self.angle_mean_ah = 20
        self.cvl_mean_bs = 1.4
        self.cvl_mean_ah = 2.2
        self.thresh_scores = 1.5

    @staticmethod
    def __scoring_tick_fn(u, mean, topx, v=0.9, p=0.5, n=0.5):
        r = np.abs((u - mean) / mean) ** n if u <= mean else (((((u - mean) * ((mean * v))) / (
        mean + ((topx - mean) * p) - mean)) / mean)) ** n if u <= mean + ((topx - mean) * p) else (v + (((((u - (
        mean + (topx - mean) * p)) * (mean - (mean * v))) / (topx - (mean + (topx - mean) * p) + (
        mean * v))) / mean))) ** n
        return r

    def _check_unified_score(self, cv1, cv2, validate=("bs", "co"), min_num_bs=1.0, dizio3d=None):

        if cv1.res2 != cv2.res1:
            return False

        sup = (cv1.res2.resseq, cv2.res2.resseq)

        cvs = (cv1, cv2), (cv2, cv1)

        for CV1, CV2 in cvs:
            if CV1.ss1 in validate or dizio3d is None:
                beta_score = self.__scoring_tick_fn(CV1.d, self.cvl_mean_bs, 2.4, v=0.9, p=0.8, n=1) + self.__scoring_tick_fn(CV2.d, self.angle_mean_bs, 180.0, v=0.9, p=0.5, n=0.5)
                alpha_score = self.__scoring_tick_fn(CV1.d, self.cvl_mean_ah, 2.4, v=0.0, p=1.0, n=1) + self.__scoring_tick_fn(CV2.d, self.angle_mean_ah, 180.0, v=0.9, p=0.5, n=0.5)


                if dizio3d is not None:
                    pass

                CV1._unified_score = np.abs(alpha_score - beta_score)

                if beta_score_uno < alpha_score and beta_score <= self.thresh_scores and CV1._unified_score >= self.params["strictness_bs"]:
                    CV1.ss2 = "bs"
                elif dizio3d is None and alpha_score < beta_score and alpha_score <= self.thresh_scores and CV1._unified_score >= self.params["strictness_ah"]:
                    CV1.ss2 = "ah"
                elif dizio3d is not None and alpha_score < beta_score and alpha_score <= self.thresh_scores and CV1._unified_score >= self.params["strictness_ah"]:
                    CV1.ss2 = "co"
                else:
                    CV1.ss2 = "co"
                if CV1.ss1 != CV1.ss2:
                    print(CV1, CV1.ss1, "-->", CV1.ss2)

        return True





class ALEPH2(_ALEPH):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    def _path_to_entity(self, entity_or_path):

        if not isinstance(entity_or_path, BIEntity):
            e = BIEntity.from_file(entity_or_path)
        else:
            e = entity_or_path
        return e

    def generate_fragments(self, entity_or_path):
        entity = self._path_to_entity(entity_or_path)
        entity = self.calculate_secondary_structure(entity)


        return entity

    def calculate_secondary_structure(self, entity_or_path):
        entity = self._path_to_entity(entity_or_path)
        log(1, "Annotating with ALEPH2:", entity)
        entity.cvectors()

        cv1 = None
        for cv2 in entity.cvectors():
            if cv1 is None:
                cv1 = cv2
                continue
            #print(cv1, cv2)
            self._check_unified_score(cv1, cv2)
            cv1=cv2


        return entity



