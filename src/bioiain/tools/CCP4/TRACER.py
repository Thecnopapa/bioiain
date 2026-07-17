import os, sys, json

import subprocess
from ...utilities.logging import log
from ...utilities.exceptions import *
from ...utilities import *

from . import CCP4_PATH





class TRACER(object):
    def __init__(self,  name=None, tracer_command="tracer"):
        self.name = name
        if name is None:
            name = ""
        else:
            name = f"_{name}"
        self.command = tracer_command
        os.makedirs(os.path.join(TEMP_FOLDER, "tracer"), exist_ok=True)
        self.instruction_file = os.path.join(TEMP_FOLDER, "tracer", f"instructions{name}.txt")



    def _write_instructions(self, **kwargs):
        with open(self.instruction_file, "w") as f:
            f.write(f"TITLE {self.name}\n")
            if kwargs.get("d", None) is not None:
                f.write(f"DEL {kwargs.get('d'):2.1f}\n")

            if kwargs.get("cell", None) is not None:
                f.write(f"CELL {' '.join([f'{c:.2f}' for c in kwargs.get('cell')])}\n")

            if kwargs.get("rcell", None) is not None:
                f.write(f"RCELL {' '.join([f'{c:.2f}' for c in kwargs.get('rcell')])}\n")

            if kwargs.get("centered", None) is not None:
                f.write(f"{kwargs.get('centered')}CENTERED\n")

            if kwargs.get("reduce"):
                f.write("REDUCED\n")
            else:
                f.write("CALCULATE\n")
        log(2, "Tracer instructions writtent to:", self.instruction_file)

    def _run_instructions(self):
        pass




    def run(self, cell:list[float], d:float=1.0,  centered:str|None=None, rcell:list|None=None, reduce:bool=True):

        self._write_instructions(cell=cell, d=d, centered=centered, rcell=rcell, reduce=reduce)
        self._run_instructions()
        return self













