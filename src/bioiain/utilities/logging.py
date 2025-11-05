
import os, shutil
import warnings


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
            elif v > 0:
                if level == 0 or level is None:
                    print(*args, **kwargs)
                elif level == "title":
                    tprint(*args, **kwargs)
                elif level == "header":
                    sprint(*args, **kwargs)
                elif type(level) is int:
                    if v >= level:
                        print1(*args, space=2*level, **kwargs)
                else:
                    print("Unknown log level: {}".format(level))



def tprint(*strings:str, head:int=10, style:str="#", end:str="\n", sep:str=" "):  # Print section title
    width = shutil.get_terminal_size()[0] -2
    string = " ".join(strings)
    tail = width - head - len(string)
    out = "\n{}{}{}{}{}".format(style*head, sep, string, sep, style*tail)
    print(out, end=end)

def eprint(*strings, head=10, style = "^", sep=" "):  # Print end of section
    tprint(*strings, head=head, style=style, end="\n\n", sep=sep)

def sprint(*strings:str, **kwargs): # Print Subtitle
    str_strings = map(str, strings)
    prefix = "\n"
    out = " # "+ " ".join(str_strings)
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






