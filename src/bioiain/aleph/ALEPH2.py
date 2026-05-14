import os, sys, json
from ..base import *
from ..utilities.exceptions import *
from ..utilities.logging import log


class ALEPH2(object):
    def __init__(self, params={}):
        self.params = dict(
            min_ah=4, min_bs=3,
            strictness_ah=0.45, strictness_bs=0.20,
            weight = "distance_avg",
        )
        self.params |= params


    def annotate(self, entity_or_path):

        if not isinstance(entity_or_path, BIEntity):
            entity = BIEntity.from_file(entity_or_path)
        else:
            entity = entity_or_path

        log(1, "Annotating with ALEPH2:", entity)

        entity.cvectors()

        pass


    def _generate_matrix_and_graph(self, entity):
        pass
