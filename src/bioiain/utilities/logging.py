
import os, sys, shutil
import warnings
import time, datetime
import math


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
    if type(level) is str:
        level = level.lower()
    if v > -2:
        if level == "error":
            if isinstance(kwargs.get("error", None), Exception):
                raise kwargs.get("error")

            elif kwargs.get("raise_exception", False):
                raise Exception(" ".join(args))
            else:
                print("\033[91m")
                print("ERROR: ", end="")
                print(*args, **kwargs)
                print("\033[0m")

        elif v > -1:
            if level == "warning":
                print("\033[93m")
                print("WARNING: ", end="")
                print(*args, **kwargs)
                print("\033[0m")
            elif level == "debug":
                print(*args, **kwargs)
            elif level == "title":
                print("\033]0;")
                print(*args, **kwargs)
                print("\a")
            elif v > 0:
                if level == 0 or level is None:
                    print(*args, **kwargs)
                elif level == "start":
                    tprint(*args, **kwargs)
                elif level == "header":
                    sprint(*args, **kwargs)
                elif level == "end":
                    eprint(*args, **kwargs)
                elif type(level) is int:
                    if v >= level:
                        print1(*args, space=2*level, **kwargs)
                    else:
                        print1("...", space=2 * level, **kwargs)
                else:
                    print("Unknown log level: {}".format(repr(level)))
                    print(*args, **kwargs)

start_time = None

def tprint(*strings:str, head:int=10, style:str="#", end:str="\n", sep:str=" ", reset_timer=True, print_timer=False):  # Print section title
    global start_time
    width = shutil.get_terminal_size()[0] -2
    string = " ".join(strings)
    timer = ""
    if print_timer and start_time is not None:
        timer = "{}{}{}".format(sep, datetime.timedelta(seconds =time.time() - start_time), sep)

    tail2 = 3 * style
    tail1_len = width - head - len(string) - len(timer) - len(tail2)

    if tail1_len < 0:
        tail1_len = 0
        tail2=""
    tail1 = style*tail1_len

    out = "\n{}{}{}{}{}{}{}".format(style*head, sep, string, sep, tail1, timer, tail2 )
    print(out, end=end)
    if reset_timer:
        start_time = time.time()

def eprint(*strings, style = "^", **kwargs):  # Print end of section
    tprint(*strings, style=style, end="\n\n", print_timer=True, **kwargs)




def sprint(*strings:str, **kwargs): # Print Subtitle
    str_strings = map(str, strings)
    prefix = "\n"
    out = " * "+ " ".join(str_strings)
    print(prefix+out,**kwargs)

def print1(*strings:str, space:int=2, **kwargs): # Print with 1 indent
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

