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

    def _check_unified_score(self, cv1, cv2, validate=("bs", "co"), min_num_bs=1.0, dizio3d=None , check_second=False):
        log(2,"Checking unified score...")

        if cv1.res2 != cv2.res1:
            return False

        sup = (cv1.res2.resseq, cv2.res2.resseq)

        if check_second:
            cvs = (cv1, cv2), (cv2, cv1)
        else:
            cvs = (cv1, cv2),

        for CV1, CV2 in cvs:
            if CV1.ss1 in validate or dizio3d is None:
                beta_score = self.__scoring_tick_fn(CV1.d, self.cvl_mean_bs, 2.4, v=0.9, p=0.8, n=1) + self.__scoring_tick_fn(CV2.d, self.angle_mean_bs, 180.0, v=0.9, p=0.5, n=0.5)
                alpha_score = self.__scoring_tick_fn(CV1.d, self.cvl_mean_ah, 2.4, v=0.0, p=1.0, n=1) + self.__scoring_tick_fn(CV2.d, self.angle_mean_ah, 180.0, v=0.9, p=0.5, n=0.5)


                if dizio3d is not None:
                    raise Exception("dizio3d is not implemented")
                    if uno[3] not in dizio3d or len(dizio3d[uno[3]]) == 0:
                        listuno = []
                        beta_score_uno += 5
                    else:
                        listuno = sorted(dizio3d[uno[3]], key=lambda x: x[2])
                        beta_score_uno += __scoring_tick_fn(listuno[0][2], dist_mean_bs, 10.0, v=0.9, p=0.5, n=0.5) if \
                        listuno[0][2] >= dist_mean_bs else 0

                    if len(listuno) < min_num_bs:
                        beta_score_uno += 5
                    else:
                        listuno = sorted(listuno, key=lambda x: x[2], reverse=True)
                        listunobleah = sorted([p[2] for p in listuno], reverse=True)
                        f = [((-1.0 * (t + 1)) / (e * len(listuno))) * 2 * (
                                    max(BS_UD_EA[int(round(listuno[t][1]))], BS_UU_EA[int(round(listuno[t][1]))]) /
                                    BS_MAX[numpy.argmax(
                                        [BS_UD_EA[int(round(listuno[t][1]))], BS_UU_EA[int(round(listuno[t][1]))]])])
                             for t, e in enumerate(listunobleah)]

                        beta_score_uno += sum(f)
                    beta_score_uno /= 4
                    alpha_score_uno /= 2

                CV1._unified_score = np.abs(alpha_score - beta_score)
                old_ss2 = CV1.ss2

                if beta_score < alpha_score and beta_score <= self.thresh_scores and CV1._unified_score >= self.params["strictness_bs"]:
                    CV1.ss2 = "bs"
                elif dizio3d is None and alpha_score < beta_score and alpha_score <= self.thresh_scores and CV1._unified_score >= self.params["strictness_ah"]:
                    CV1.ss2 = "ah"
                elif dizio3d is not None and alpha_score < beta_score and alpha_score <= self.thresh_scores and CV1._unified_score >= self.params["strictness_ah"]:
                    CV1.ss2 = "co"
                else:
                    CV1.ss2 = "co"
                if old_ss2 is not None:
                    log(3,CV1, CV1.ss1, "-->", old_ss2, "-->", CV1.ss2)
                    if old_ss2 != CV1.ss2:
                        raise Exception("ss2 calculations do not match")
                if CV1.ss1 != CV1.ss2:
                    log(3, CV1, CV1.ss1, "-->", CV1.ss2)

        return True

    def _check_impossible_angle(self, cv, cv3, cv2, cv1):
        log(2, "Checking impossible angle...")
        print(cv, cv3)
        print(cv.ss2, cv3.ss2)
        if cv.res1.resnum -1 != cv3.res3.resnum:
            log("warning", f"CVs not continuous ({cv3.res3.resnum}-{cv.res1.resnum})")
            return False

        impossible = False
        # TODO: check this
        if cv.ss2 == "bs":
            if abs(cv3.a - self.angle_mean_bs) > 50:
                impossible = True
        elif cv.ss2 == "ah":
            if abs(cv3.a - self.angle_mean_ah) > 50:
                impossible = True
        if impossible:
            pass










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

        cvectors = entity.cvectors()
        nvectors = len(cvectors)
        cv1 = None
        for n, cv2 in enumerate(cvectors):
            if cv1 is None:
                cv1 = cv2
                continue
            self._check_unified_score(cv1, cv2, check_second=n==nvectors-1)
            if n > 3:
                self._check_impossible_angle(cv1, cvectors[n-4])
            cv1=cv2


        return entity
