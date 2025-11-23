import os, sys, subprocess
from idlelib.undo import Command

from ..utilities import *
from ..biopython.base import BiopythonOverlayClass



pymol_colours = ['green', 'cyan', 'red', 'yellow', 'violet','blue',
               'salmon', 'lime', 'pink', 'slate', 'magenta', 'orange', 'marine',
               'olive', 'purple', 'teal', 'forest', 'firebrick', 'chocolate',
               'wheat', 'white', 'grey']



def quick_display(entity:BiopythonOverlayClass|list[BiopythonOverlayClass]) -> str:
    """
    Displays entity or list of entities with PyMol. Exports entities to ./.temp and saves generated script in the same
    directory as quick_display.pml . Entities are named as N_[entity_id] following input order.
    :param entity: Entity or list of entities.
    :return: Path to the generated script.
    """
    script = PymolScript("quick_display")
    os.makedirs("./.temp", exist_ok=True)

    if type(entity) is not list:
        entity = [entity]
    for n, entity in enumerate(entity):
        name = "{}_{}".format(n, entity.id)
        script.load_entity(entity, name, overwrite=False)
    script.write_script("./.temp")
    script.execute()
    return script.path






class PymolScript(object):
    """
    Class to build PyMol scripts from predetermined functions or custom ones.
    :param name: Name of the script. Will de set as a filename. Default is ".temp_pymol_script".
    :return: PymolScript Object.
    """
    def __init__(self, name="temp_pymol_script", folder:str="./temp", pymol_path = "pymol"):
        self.pymol_path = pymol_path
        self._bioiain = "bioiain"
        self.name = name
        self.folder = folder
        os.makedirs(self.folder, exist_ok=True)
        self.subfolder = os.path.join(folder, self.name)
        os.makedirs(self.subfolder, exist_ok=True)
        self.input = {}
        self.data = {}
        self.commands = []
        self.path = None


    def write_script(self, filename:str=None) -> str:
        """
        Writes the stored commands to a file. The file can be executed from the terminal or run as a PyMol script.
        :param filename: (optional) Path to the file to write to. Uses script name and current wd as default.
        :return: Path to the file.
        """
        if filename is None:
            filename = self.name
        filepath = os.path.join(self.subfolder, filename+".pml")
        with open(filepath, "w") as f:
            f.write("pymol\n\n")
            #f.write("import {} as bi\n\n\n".format(self._bioiain))
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
            self.write_script(".")
        cmd = [self.pymol_path, self.path]

        logging.log("debug", "$ " + " ".join(cmd))
        subprocess.run(cmd)


    def add(self, fun, *args, **kwargs) -> Command:
        """
        Generates a command object and adds the command to the script. Can be used to insert custom functions. Beware when adding parameters as strings.
        must be double-quoted ("'string'").
        :param fun: Name(string) of function to execute. Use .import() to import such function if necessary.
        :param args: Args to pass the function. "strings" -> variables, "'strings'" -> strings.
        :param kwargs: Same as args but with keywords.
        :return: Generated Command object.
        """
        c = self.Command(fun, *args, **kwargs)
        self.commands.append(c)
        return c



    class Command(object):
        """
        Class for commands stored PyMol script.
        :param fun: Name(string) of function to execute. Use .import() to import such function if necessary.
        :param args: Args to pass the function. "strings" -> variables, "'strings'" -> strings.
        :param to: Name of variable to assign the return of the function.
        :param is_cmd: whether the command is within od pymol.cmd.
        :param kwargs: Same as args but with keywords.
        :return: Command object.
        """
        def __init__(self, fun:str, *args, to:str=None, is_cmd=True, **kwargs):
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


        def construct_command(self) -> str:
            """
            Generates final string to append to script. Uses parameters stored in the instance.
            :return: Generated string.
            """
            if type(self.args) == str:
                arg_str = self.args
            else:
                arg_str = ", ".join(self.args)
            kwarg_str = ", ".join([f"{k}={v}" for k, v in self.kwargs.items()])
            c = "{}({},{})".format(self.fun, arg_str, kwarg_str)
            if self.is_cmd:
                c = "cmd."+ c
            if self.to is not None:
                c = "{} = {}".format(self.to, c)
            self.cmd = c
            return self.cmd

    @staticmethod
    def _process_sele(sele:str, force_str:bool=False) -> str:
        """
        Adds extra quotes to PyMol selections. Selections are strings within brackets e.g (all).
        :param sele: Selection string.
        :param force_str: Whether to always double-quote input string.
        :return: Double-quoted selection or string if not a selection.
        """
        if sele.startswith("(") and sele.endswith(")") or force_str:
            sele = f"'{sele}'"
        else:
            sele = f"{sele}"
        return sele

    @staticmethod
    def _to_str(string):
        return f"'{string}'"


    def print(self, *args, **kwargs) -> Command:
        """
        Adds command to print with builtin print.
        :param args: Args to pass to print.
        :param kwargs: Kwargs to pass to print.
        :return: Generated Command object. -> Nothing
        """
        fun = "print"
        return self.add(fun, *args, is_cmd=False, **kwargs)


    def load(self, path:str, name:str, **kwargs) -> Command:
        """
        Adds command to load file from path.
        :param path: Path to file.
        :param name: Name of created PyMol object.
        :param kwargs:
        :return: Generated Command object -> Unknown.
        """
        fun = "load"
        args = f"'{os.path.abspath(path)}'", f"'{name}'"
        return self.add(fun, *args, **kwargs)

    def load_entity(self, entity:BiopythonOverlayClass, name:str|None=None, overwrite:bool=True) -> Command:
        """
        Adds command to load file from entity. Entity is exported to ./.temp as of the cwd.
        :param entity:
        :param name: (optional) Name of created PyMol object. Defaults to the entity id
        :return: Generated Command object -> Unknown.
        """
        if name is None:
            name = entity.data["info"]["name"]
        if not overwrite:
            n = 1
            while name+".pdb" in os.listdir(self.subfolder):
                name = "{}_{}".format(entity.data["info"]["name"], n)
                n += 1
        path = entity.export(self.subfolder, name, data=True)[0]
        return self.load(path, name)


    def disable(self, sele:str, **kwargs) -> Command:
        """
        Adds Command to disable selection.
        :param sele: Selection string.
        :param kwargs:
        :return: Generated Command object -> Unknown.
        """
        sele = self._process_sele(sele, **kwargs)
        fun = "disable"
        return self.add(fun, sele, **kwargs)


    def symmetries(self, obj:str="original", prefix:str="sym", distance:int=6, **kwargs) -> Command:
        fun = "symexp"
        obj = self._to_str(obj)
        args = [self._to_str(prefix), obj, obj, str(distance)]
        return self.add(fun, *args, **kwargs)


    def cell(self, **kwargs) -> Command:
        fun = "show"
        args = "'cell'"
        return self.add(fun, args, **kwargs)

    def group(self, prefix:str="sym", name:str|None=None, **kwargs) -> Command:
        fun = "group"
        sele = self._to_str(prefix+"*")
        if name is None:
            name = prefix
        args = [self._to_str(name), sele]
        return self.add(fun, *args, **kwargs)


    def color(self, sele:str, color:str|int="black", **kwargs) -> Command:
        fun = "color"
        sele = self._to_str(sele)
        if type(color) is int:
            color = pymol_colours[len(pymol_colours) % (color+1)]
        color = self._to_str(color)
        args = [color, sele]

        return self.add(fun, *args, **kwargs)

    def spectrum(self, sele: str, spectrum: str = "b", color="rainbow", **kwargs) -> Command:
        fun = "spectrum"
        sele = self._to_str(sele)
        spectrum = self._to_str(spectrum)
        color = self._to_str(color)
        args = [spectrum, color, sele]
        #kwargs["spectrum"] = spectrum

        return self.add(fun, *args, **kwargs)





    """def pymol_group(identifier="sym", name=None, quiet=False):
        group = []
        if name is None:
            name = identifier
        if not quiet:
            print("(PyMol) Grouping:", identifier, "in", name)
        for obj in pymol.cmd.get_names(type='objects'):
            if type(identifier) is str:
                if identifier in obj:
                    group.append(obj)
            if type(identifier) is list:
                for i in identifier:
                    if i in obj:
                        group.append(obj)

        pymol.cmd.group(name, " ".join(group))"""













