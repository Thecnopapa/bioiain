import os, sys, subprocess
from ..utilities import *

class PymolScript(object):
    def __init__(self, name=".temp_pymol_script", pymol_path = "pymol"):
        self.pymol_path = pymol_path
        self._bioiain = "bioiain"
        self.name = name
        self.input = {}
        self.data = {}
        self.commands = []
        self.path = None


    def write_script(self, filepath=None):
        if filepath is None:
            filepath = self.name+".py"
        with open(filepath, "w") as f:
            f.write("#!pymol\n\n")
            f.write("import {} as bi\n\n\n".format(self._bioiain))
            for cmd in self.commands:
                f.write(repr(cmd)+"\n")
        self.path = os.path.abspath(filepath)


    def execute(self):
        if self.path is None:
            self.write_script()
        subprocess.run([self.pymol_path, self.path])


    def add(self, fun, *args, **kwargs):
        c = self.Command(fun, *args, **kwargs)
        self.commands.append(c)
        return c



    class Command(object):
        def __init__(self, fun, *args, to=None, is_cmd=True, **kwargs):
            self.fun = fun
            self.args = args
            self.kwargs = kwargs
            self.to = to
            self.is_cmd = is_cmd
            self.cmd = None

        def __repr__(self):
            if self.cmd is None:
                self.construct_command()
            return self.cmd


        def construct_command(self):
            if type(self.args) == str:
                arg_str = self.args
            else:
                arg_str = ", ".join(self.args)
            kwarg_str = ", ".join(["{}={}" for k, v in self.kwargs.items()])
            c = "{}({},{})".format(self.fun, arg_str, kwarg_str)
            if self.is_cmd:
                c = "cmd."+ c
            if self.to is not None:
                c = "{} = {}".format(self.to, c)
            self.cmd = c
            return self.cmd

    def _process_sele(self, sele, force_str=False, **kwargs):
        if sele.startswith("(") and sele.endswith(")") or force_str:
            sele = f"'{sele}'"
        else:
            sele = f"{sele}"
        return sele


    def print(self, *args, **kwargs):
        fun = "print"
        return self.add(fun, *args, is_cmd=False, **kwargs)


    def load(self, path:str, name:str, **kwargs):
        fun = "load"
        args = f"'{path}'", f"'{name}'"
        return self.add(fun, *args, **kwargs)


    def disable(self, sele:str, **kwargs):
        sele = self._process_sele(sele, **kwargs)
        fun = "disable"
        return self.add(fun, sele, **kwargs)
















