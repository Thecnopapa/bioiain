import os, sys, subprocess, shutil

from ..base.atom import PseudoAtom
from ..utilities import relative_path
from ..utilities.logging import log
from ..utilities import *



pymol_colours = ['green', 'cyan', 'red', 'yellow', 'violet','blue',
               'salmon', 'lime', 'pink', 'slate', 'magenta', 'orange', 'marine',
               'olive', 'purple', 'teal', 'forest', 'firebrick', 'chocolate',
               'wheat', 'white', 'grey']



def quick_display(entity, execute=True) -> str:
    """
    Displays entity or list of entities with PyMol. Exports entities to ./.temp and saves generated script in the same
    directory as quick_display.pml . Entities are named as N_[entity_id] following input order.
    :param entity: Entity or list of entities.
    :return: Path to the generated script.
    """
    script = PymolScript("quick_display", folder=os.path.join(TEMP_FOLDER,"pml"))

    if type(entity) is not list:
        entity = [entity]
    for n, entity in enumerate(entity):
        name = "{}_{}".format(n, entity.name())
        script.load_entity(entity, name, overwrite=False)
    script.write_script()
    if execute:
        script.execute()
    return script






class PymolScript(object):
    """
    Class to build PyMol scripts from predetermined functions or custom ones.
    :param name: Name of the script. Will de set as a filename. Default is ".temp_pymol_script".
    :return: PymolScript Object.
    """
    def __init__(self, name="pymol_script", folder:str|None=None, tmp_folder=None, pymol_path = "pymol", use_temp=False, allow_plugins=False):
        self.pymol_path = pymol_path
        self._bioiain = "bioiain"
        self.name = name
        if folder is None:
            folder = os.path.join(SUBDIR_NAME, "pml_sessions")
        self.folder = folder
        if tmp_folder is None:
            tmp_folder = os.path.join(TEMP_FOLDER, "pml")
        self.tmp_folder=tmp_folder
        os.makedirs(self.folder, exist_ok=True)
        os.makedirs(self.tmp_folder, exist_ok=True)
        self.use_temp = use_temp
        if use_temp:
            self.subfolder =  os.path.join(self.tmp_folder, self.name)
        else:
            self.subfolder = os.path.join(self.folder, self.name)

        os.makedirs(self.subfolder, exist_ok=True)
        self.input = {}
        self.data = {}
        self.commands = []
        self.path = None
        self.session_path = None
        log(1, f"New: {self}")


    def __repr__(self):
        if self.path is None:
            return f"<bi.{self.__class__.__name__}: {self.name} UNSAVED in folder: {self.subfolder} N={len(self.commands)}>"
        else:
            return f"<bi.{self.__class__.__name__}: {self.name} N={len(self.commands)}>"


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
        def __init__(self, fun:str, *args, to:str=None, is_cmd=True, is_fun=True, **kwargs):
            self.fun = fun
            self.args = args
            self.kwargs = kwargs
            self.to = to
            self.is_cmd = is_cmd
            self.is_fun = is_fun
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
            arglist = []
            if len(arg_str) > 0:
                arglist.append(arg_str)
            if len(kwarg_str) > 0:
                arglist.append(kwarg_str)
            if self.is_fun:
                c = "{}({})".format(self.fun, ", ".join(arglist))
            else:
                c = "{} {}".format(self.fun, ", ".join(arglist))
            if self.is_cmd:
                c = "cmd."+ c
            if self.to is not None:
                c = "{} = {}".format(self.to, c)
            self.cmd = c
            return self.cmd

    def write_script(self, *args, **kwargs) -> str:
        return self.save(*args, **kwargs)

    def save(self, filename:str=None) -> str:
        """
        Writes the stored commands to a file. The file can be executed from the terminal or run as a PyMol script.
        :param filename: (optional) Path to the file to write to. Uses script name and current wd as default.
        :return: Path to the file.
        """
        if filename is None:
            filename = self.name
        filepath = os.path.join(self.subfolder, filename+".script.pml")
        with open(filepath, "w") as f:
            f.write(f"#!{self.pymol_path}\n\n")
            #f.write(f"try:\n\tcd os.path.dirname(sys.argv[1])\nexcept:\n\tpass\n")
            #f.write("import os\n")
            f.write(f"os.chdir('{self.subfolder}')\n")

            #f.write("import {} as bi\n\n\n".format(self._bioiain))
            for cmd in self.commands:
                f.write(repr(cmd)+"\n")
        self.path = os.path.abspath(filepath)
        log(1, f"PyMol Session saved at:\npymol {self.path}")
        try:
            os.chmod(self.path, 0o755)
        except:
            log("warning", "Could not give exec permissions to script")
        return self.path

    def compile(self, **kwargs):
        return self.execute(compile=True, **kwargs)

    def execute(self, quiet=True, pymol_path=None, full_screen=False, compile=False, use_compiled=True, extra_options="-k"):
        """
        Executes the script on the current thread. Not sure if it is blocking or not.
        """
        if self.path is None:
            self.write_script()
        if pymol_path is None:
            pymol_path = self.pymol_path
        cmd = [pymol_path]
        if extra_options is not None:
            if type(extra_options) is str:
                if len(extra_options) > 0:
                    cmd.append(extra_options)
            elif type(extra_options) is list:
                cmd.extend(extra_options)

        if full_screen:
            cmd.extend(["-x", "-e"])

        if quiet:
            cmd.extend(["-qQ"])

        tmp_session_path = None
        if compile:
            self.session_path = self.path.replace(".script.pml", ".session.pse")
            tmp_session_path = self.path.replace(".script.pml", ".compiler.pml")
            shutil.copy(self.path, tmp_session_path)
            with open(tmp_session_path, "a") as f:
                f.write(f"\ncmd.save('{self.session_path}')")
            cmd.extend(["-c"])
            cmd.extend(["-l", tmp_session_path])
        else:
            if use_compiled and self.session_path is not None:
                cmd.extend([self.session_path])
            else:
                cmd.extend(["-l", self.path])

        logging.log("debug", "$ " + " ".join(cmd))
        try:
            subprocess.run(cmd, cwd=self.subfolder)
            if compile:
                log(1, f"PyMol Session compiled at: pymol {self.session_path}")
                os.remove(tmp_session_path)
        except KeyboardInterrupt:
            logging.log("debug", "\nClosing Pymol...")
        except Exception as e:
            log("error","(PYMOL)", e)


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

    def raw(self, line, *args, **kwargs):
        if type(line) is str:
            line = [line]
        self.add(", ".join(line), is_cmd=False, is_fun=False *args, **kwargs)
        return self





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
            sele = f"'{sele}'"
        return sele

    @staticmethod
    def _to_str(string) -> str:
        return f"'{string}'"


    def print(self, *args, literal=True, **kwargs) -> str:
        """
        Adds command to print with builtin print.
        :param args: Args to pass to print.
        :param kwargs: Kwargs to pass to print.
        :return: Generated Command object. -> Nothing
        """
        fun = "print"
        if literal:
            new_args = []
            for arg in args:
                new_args.append(self._to_str(repr(arg).replace("'", "\"")))
            args = new_args
        self.add(fun, *args, is_cmd=False, **kwargs)
        return " ".join(args)


    def load(self, path:str, name:str=None, create=False, **kwargs) -> str:
        """
        Adds command to load file from path.
        :param path: Path to file.
        :param name: Name of created PyMol object.
        :param kwargs:
        :return: Generated Command object -> Unknown.
        """
        if create:
            fun = "create"
        else:
            fun = "load"

        path = relative_path(path, self.subfolder)
        if name is None:
            name = os.path.basename(path).split(".")[0]

        args = f"'{path}'", f"'{name}'"
        self.add(fun, *args, **kwargs)
        return name

    def delete(self, sele="(all)", **kwargs):

        fun = "delete"
        sele = self._to_str(sele)
        self.add(fun, sele, **kwargs)
        return self


    def orient(self, sele="(all)", **kwargs):

        fun = "orient"
        sele = self._to_str(sele)
        self.add(fun, sele, **kwargs)
        return self

    def center(self, sele="(all)", **kwargs):

        fun = "center"
        sele = self._to_str(sele)
        self.add(fun, sele, **kwargs)
        return self

    def show(self, sele="(all)", representation="cartoon", **kwargs):
        fun = "show"
        sele = self._to_str(sele)
        self.add(fun, self._to_str(representation), sele, **kwargs)
        return self

    def hide(self, sele="(all)", representation="everything", **kwargs):
        fun = "hide"
        sele = self._to_str(sele)
        self.add(fun, self._to_str(representation), sele, **kwargs)
        return self

    def load_entity(self, entity, name:str|None=None, overwrite:bool=True, minimal=True) -> str:
        """
        Adds command to load file from entity. Entity is exported to t/mp/bioiain/pymol as of the cwd.
        :param entity:
        :param name: (optional) Name of created PyMol object. Defaults to the entity id
        :return: Generated Command object -> Unknown.
        """
        if name is None:
            name = entity.data["info"]["name"]
        if not overwrite:
            n = 1
            while name+".cif" in os.listdir(self.subfolder):
                name = "{}_{}".format(entity.data["info"]["name"], n)
                n += 1

        folder = self.subfolder
        path = entity.export(minimal=minimal, target_folder=folder)
        log(2, "Entity saved to:", path)
        self.load(path, name)
        return name


    def disable(self, sele:str, **kwargs):
        """
        Adds Command to disable selection.
        :param sele: Selection string.
        :param kwargs:
        :return: Generated Command object -> Unknown.
        """
        sele = self._to_str(sele)
        fun = "disable"
        self.add(fun, sele, **kwargs)
        return self


    def symmetries(self, obj:str="original", prefix:str="sym", distance:int=6, **kwargs):
        fun = "symexp"
        obj = self._to_str(obj)
        args = [self._to_str(prefix), obj, obj, str(distance)]
        self.add(fun, *args, **kwargs)
        return self


    def cell(self, **kwargs):
        fun = "show"
        args = "'cell'"
        self.add(fun, args, **kwargs)
        return self

    def group(self, prefix:str="sym", name:str|None=None, also_suffix=False, **kwargs) -> str:
        fun = "group"
        if also_suffix:
            prefix = "*"+prefix if not prefix.startswith("*") else prefix
        sele = self._to_str(prefix+"*" if  not prefix.endswith("*") else prefix)
        if name is None:
            name = prefix
        name = self._to_str(name)
        args = [name, sele]
        self.add(fun, *args, **kwargs)
        return name

    def align(self, moving:str, fixed:str, fun="align", **kwargs):
        sele_fixed = self._to_str(fixed)
        sele_moving = self._to_str(moving)

        args = [sele_moving, sele_fixed]
        self.add(fun, *args, **kwargs)
        return self

    def merge(self, target:str, sele:str, state=-1, **kwargs) -> str:
        fun = "create"
        sele_target= self._to_str(target)
        sele = self._to_str(sele)

        args = [sele_target, sele]
        self.add(fun, *args, state=state, **kwargs)
        return sele_target


    def color(self, sele:str, color:str|int="black", **kwargs):
        fun = "color"
        sele = self._to_str(sele)
        if type(color) is int:
            color = pymol_colours[len(pymol_colours) % (color+1)]
        color = self._to_str(color)
        args = [color, sele]

        self.add(fun, *args, **kwargs)
        return self

    def spectrum(self, sele:str="(all)", spectrum: str = "b", color="rainbow", **kwargs):
        fun = "spectrum"
        sele = self._to_str(sele)
        spectrum = self._to_str(spectrum)
        color = self._to_str(color)
        args = [spectrum, color, sele]
        #kwargs["spectrum"] = spectrum

        self.add(fun, *args, **kwargs)
        return self

    def pseudoatom(self, name="tmp", coord=(0,0,0), atomname="PA", **kwargs) -> str:
        fun = "pseudoatom"
        name = self._to_str(name)
        self.add(fun, name, pos=coord, name=atomname **kwargs)
        return name

    def line(self, name="line", sele1=None, sele2=None, coord1=(0,0,0), coord2=(0,0,0), show_distance=False, **kwargs) -> str:
        fun = "distance"
        if isinstance(coord1, PseudoAtom):
            coord1 = coord1.coord
        if isinstance(coord2, PseudoAtom):
            coord2 = coord2.coord
        if sele1 is None:
            sele1 = "tmp1"
            self.pseudoatom(sele1, coord=[float(c) for c in coord1], **kwargs)
        if sele2 is None:
            sele2 = "tmp2"
            self.pseudoatom(sele2, coord=[float(c) for c in coord2], **kwargs)
        args = [self._to_str(name), self._to_str(sele1), self._to_str(sele2)]

        r = self.add(fun, *args, **kwargs)
        if sele1 == "tmp1":
            self.delete(sele1)
        if sele2 == "tmp2":
            self.delete(sele2)
        if not show_distance:
            self.hide(name, "label")
        return name

    def set(self, param, value):
        fun = "set"
        param = self._to_str(param)
        value = self._to_str(value)
        self.add(fun, param, value)
        return self

















