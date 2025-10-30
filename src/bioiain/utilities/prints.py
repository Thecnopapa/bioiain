
import sys, os, shutil




def print(*args, **kwargs):
    try:
        if os.environ["BI_VERBOSE"]:
            vars.log(*args, **kwargs)
    except:
        pass
    import builtins
    builtins.print(*args, **kwargs)


def tprint(*strings, head=10, style="#", end="\n", sep=" ", log=True):  # Print section title
    width = shutil.get_terminal_size()[0] -2
    string = " ".join(strings)
    tail = width - head - len(string)
    out = "\n{}{}{}{}{}".format(style*head, sep, string, sep, style*tail)
    try:
        if not vars.verbose:
            vars.log(out, end = end, timestamp=False)
    except:
        pass
    print(out, end=end)

def eprint(*strings, head=10, style = "^", sep=" "):  # Print end of section
    tprint(*strings, head=head, style=style, end="\n\n", sep=sep)

def sprint(*strings,**kwargs): # Print Subtitle
    str_strings = map(str, strings)
    prefix = "\n"
    out = " # "+ " ".join(str_strings)
    try:
        if not vars.verbose:
            vars.log(out,prefix = prefix, **kwargs)
    except:
            pass
    if globals_loaded and "quiet" in vars:
        if vars.quiet:
            return

    print(prefix+out,**kwargs)

def print1(*strings, space=2, log=True, **kwargs): # Print with 1 indent
    str_strings = []
    for string in strings:
        if type(string) == list or type(string) == tuple:
            for string2 in string:
                str_strings.append(str(string2))
        else:
            str_strings.append(str(string))
    #str_strings = map(str, strings)
    out = "{}> {}".format(" " * space, " ".join(str_strings))
    try:
        if not vars.verbose:
            vars.log(out, **kwargs)
    except:
        pass
    if globals_loaded:
        if vars.quiet:
            return
    print(out, **kwargs)

def print2(*strings, **kwargs): # Print with 2 indents
    print1(strings, space=4, **kwargs)

def print3(*strings, **kwargs):
    print1(strings, space=6, **kwargs)

def print4(*strings, **kwargs):
    print1(strings, space=8, **kwargs)

def print5(*strings, **kwargs):
    print1(strings, space=10, **kwargs)

def print6(*strings, **kwargs):
    print1(strings, space=12, **kwargs)