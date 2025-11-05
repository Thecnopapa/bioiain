import os, sys, subprocess
from ..utilities import *

class PymolScript(object):
    """
    Class to build PyMol scripts from predetermined functions or custom ones.
    :param name: Name of the script. Will de set as a filename. Default is ".temp_pymol_script".
    :return: PymolScript Object.
    """
    def __init__(self, name=".temp_pymol_script", pymol_path = "pymol"):
        self.pymol_path = pymol_path
        self._bioiain = "bioiain"
        self.name = name
        self.input = {}
        self.data = {}
        self.commands = []
        self.path = None


    def write_script(self, filepath:str=None) -> str:
        """
        Writes the stored commands to a file. The file can be executed from the terminal or run as a PyMol script.
        :param filepath: (optional) Path to the file to write to. Uses script name and current wd as default.
        :return: Path to the file.
        """
        if filepath is None:
            filepath = self.name+".py"
        with open(filepath, "w") as f:
            f.write("pymol\n\n")
            f.write("import {} as bi\n\n\n".format(self._bioiain))
            for cmd in self.commands:
                f.write(repr(cmd)+"\n")
        self.path = os.path.abspath(filepath)
        os.chmod(self.path, 0o755)
        return self.path


    def execute(self):
        """
        Executes the script on the current thread. Not sure if it is blocking or not.
        """
        if self.path is None:
            self.write_script()
        cmd = [self.pymol_path, self.path]

        logging.log("debug", "$ " + " ".join(cmd))
        subprocess.run(cmd)


    def add(self, fun, *args, **kwargs):
        """
        Adds command to the script. Can be used to insert custom functions. Beware when adding parameters as strings.
        must be double-quoted ("'string'").
        :param fun: Function to execute. Use .import() to import such function if necessary.
        :param args: Args to pass the function. "strings" -> variables, "'strings'" -> strings.
        :param kwargs: Same as args but with keywords.
        :return: Command object.
        """
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
















