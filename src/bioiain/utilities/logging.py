
import os, sys, shutil, time, datetime, requests

from .. import SUBDIR_NAME, TEMP_FOLDER, WD, FD

class Log(object):
    def __init__(self):
        self.stdout = None
        self.stderr = None
        self.folder = "./logs"
        self.files = {"default": os.path.join(SUBDIR_NAME,"default.log"),
                      "debug": os.path.join(SUBDIR_NAME, "debug.log")}
        self.logging = True
        self.terminal = True

    def __repr__(self):
        return f"<bi.Log: default: {self.files['default']}>"

    def list(self):
        return [l for l in self.files.values()]

    def dict(self):
        return self.files

    def __add__(self, s):
        pass

    def __call__(self):
        pass

    def log(self):
        pass

    def error(self):
        pass

    def warning(self):
        pass

    def title(self):
        pass

    def start(self):
        pass

    def end(self):
        pass

    def pause(self):
        pass

    def resume(self):
        pass

    def set_stdout(self, filepath):
        pass

    def set_stderr(self, filepath):
        pass

    def set_log_file(self, filepath, log_name="default"):
        pass

    def add_timestamp(self, log_name=None):
        pass

    def disable(self):
        self.terminal = False
        self.logging = False


def colour(colour:str|None=None, string:str|None=None, end:bool=True):
    if colour is None:
        colour="end"
    if string is None:
        return f"\033[{colour_list.get(colour, 0)}m"
    else:
        s = f"\033[{colour_list.get(colour, 0)}m{string}"
        if end:
            s+=f"\033[0m"

        return s

colour_list = {
    "end":0,
    "black":90,
    "red":91,
    "green":92,
    "yellow":93,
    "blue":94,
    "magenta":95,
    "cyan":96,
    "white":97,
}





def log(level:int|str=1, *args, **kwargs):
    """
    Log a message and display it according to the given level if higher than environment variable "BI_VERBOSE".
    Builtin prints are always displayed.
    If unset BI_VERBOSE is set to 10.
    BI_VERBOSE == 0 displays only ERROR, WARNING and DEBUG messages.
    BI_VERBOSE == -1 display only ERROR.
    BI_VERBOSE == -1 display nothing.
    :param level: Verbose level: ERROR | WARNING | DEBUG | TITLE | HEADER | int
    :param args: args for print function
    :param kwargs: kwargs for print function
    """
    v = int(os.environ.get("BI_VERBOSE", 10))
    try:
        from .dataframes import print_df
        from polars import DataFrame
        if type(level) is DataFrame:
            print_df(level, *args, **kwargs)
            return level
    except:
        pass
    if type(level) is str:
        try:
            level = int(level)
        except:
            level = level.lower()


    if v > -2:
        if level == "error":
            if isinstance(kwargs.get("error", None), Exception):
                raise kwargs.get("error")

            elif kwargs.get("raise_exception", False):
                raise Exception(" ".join([str(a) for a in args]))
            else:
                print(colour("red"))
                print("ERROR: ", end="")
                print(*args, **kwargs)
                print(colour("end"))

        elif v > -1:
            if level == "warning":
                print(colour("yellow"),end="")
                print("WARNING: ", end="")
                print(*args, **kwargs)
                print(colour("end"),end="")
            elif level == "debug":
                print(*args, **kwargs)
            elif level == "title":
                print("\033]0;",end="")
                print(*args, **kwargs)
                print("\a",end="")
            elif v > 0:
                if level == 0 or level is None:
                    print(*args, **kwargs)
                elif level == "start":
                    _tprint(*args, **kwargs)
                elif level == "header":
                    print(colour("white"), end="")
                    _sprint(*args, **kwargs)
                    print(colour("end"), end="")
                elif level == "end":
                    _eprint(*args, **kwargs)
                elif type(level) is int:
                    if v >= level:
                        _print1(*args, space=2*level, **kwargs)
                    else:
                        _print1("...", space=2 * level, **kwargs)
                else:
                    print("Unknown log level: {}".format(repr(level)))
                    print(*args, **kwargs)

start_time = None

def _tprint(*strings:str, head:int=10, style:str="#", end:str="\n", sep:str=" ", reset_timer=True, print_timer=False):  # Print section title
    global start_time
    width = shutil.get_terminal_size()[0] -2
    string = " ".join([str(s) for s in strings])
    timer = ""
    if print_timer and start_time is not None:
        timer = "{}{}{}".format(sep, datetime.timedelta(seconds =time.time() - start_time), sep)

    tail2 = 3 * style
    tail1_len = width - head - len(string) - len(timer) - len(tail2)

    if tail1_len < 0:
        tail1_len = 0
        tail2=""
    tail1 = style*tail1_len

    out = "\n{}{}{}{}{}{}{}".format(style*head, sep, colour("white", string), sep, tail1, timer, tail2 )
    print(out, end=end)
    if reset_timer:
        start_time = time.time()

def _eprint(*strings, style = "^", print_timer=True, **kwargs):  # Print end of section
    _tprint(*strings, style=style, end="\n\n", print_timer=print_timer, **kwargs)




def _sprint(*strings:str, **kwargs): # Print Subtitle
    str_strings = map(str, strings)
    prefix = "\n"
    out = " * "+ " ".join(str_strings)
    print(prefix+out,**kwargs)

def _print1(*strings:str, space:int=2, **kwargs): # Print with 1 indent
    str_strings = []
    for string in strings:
        if type(string) == list or type(string) == tuple:
            for string2 in string:
                str_strings.append(str(string2))
        else:
            str_strings.append(str(string))
    #str_strings = map(str, strings)
    out = "{}> {}".format(" " * space, " ".join(str_strings))
    print(out, **kwargs)


def print_children(d):
    if type(d) == list:
        d = d[0]
        print("(list)[0]")
    print("strings:")
    [print(k, v) for k, v in d.items() if type(v) == str]
    print("other:")
    [print(k, type(v), len(v)) for k, v in d.items() if type(v) != str and v is not None]



try:
    if original_std_out is None:
        raise Exception("")
except:
    original_std_out = sys.stdout
    std_out = sys.stdout

def change_std_out(target_file, mode="w"):
    global std_out
    sys.std_out = open(target, mode)
    std_out = sys.stdout


def restore_std_out():
    global std_out
    sys.stdout = original_std_out
    std_out = sys.stdout


def stop_logging():
    sys.std_out = open(os.devnull, "w")


def resume_logging():
    global std_out
    sys.stdout = std_out




def send_tensorboard_run(host, folder, run, file, key, epoch=0, protocol="https"):

    url = f"{protocol}://{host}/runs/"
    fname = os.path.basename(file)
    fname = fname.replace(".0", f".{epoch}")
    log("header", "Uploading run to:", url)
    assert key is not None

    with open(file, "rb") as f:

        resp = requests.put(
                url,
                headers={
                    "Content-Type": "application/x-www-form-urlencoded",
                    "key":key,
                    "folder":folder,
                    "run":run,
                    "fname":fname
                },
                data=f.read(),
                timeout=3000
                )
        print(resp.text)
        if resp.status_code != 200:
            print(resp)
            raise Exception(f"Error [{resp.status_code}] uploading file to: {url}")
    return resp



TRACE = ("--trace" in sys.argv) or ("--tracemalloc" in sys.argv)
CURRENTLY_TRACING = False

if TRACE:
    tracemalloc_start()

def tracemalloc_start():
    global CURRENTLY_TRACING
    log("header", f"Starting tracemalloc")
    import tracemalloc
    tracemalloc.start()
    CURRENTLY_TRACING = True

def tracemalloc_stop():
    global CURRENTLY_TRACING
    log("header", f"Starting tracemalloc")
    import tracemalloc
    tracemalloc.stop()
    CURRENTLY_TRACING = False

def tracemalloc_top(top=15):
    global CURRENTLY_TRACING
    if CURRENTLY_TRACING:
        import tracemalloc
        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics('lineno')

        log("header", f"Tracemalloc top {top}")
        for stat in top_stats[:top]:
            log(1, stat)
        print()
