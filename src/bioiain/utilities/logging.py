import os, sys, shutil, time, datetime, requests, curses, math
import numpy as np

from .. import SUBDIR_NAME, TEMP_FOLDER, WD, FD
CURSED = False
CURSED_LOG = None

class Log(object):
    def __init__(self):
        self.stdout = None
        self.stderr = None
        self.folder = "./logs"
        self.files = {"default": os.path.join(SUBDIR_NAME,"default.log"),
                      "debug": os.path.join(SUBDIR_NAME, "debug.log")}
        self.logging = True
        self.terminal = True
        self.cursed_log = None

    def __repr__(self):
        return f"<bi.{self.__class__.__name__}: default: {self.files['default']}>"

    def list(self):
        return [l for l in self.files.values()]

    def dict(self):
        return self.files

    def __add__(self, s):
        pass

    def __call__(self):
        pass

    def __getitem__(self, logname):
        return self.files.get(logname, None) 

    def 
    def _process_input()

    def log(self):

        with open(self["debug"]) as debug_f:
            debug_f.write()

        if not logging:
            return

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
    
    if CURSED and CURSED_LOG is not None:
        def print(*args, **kwargs):
            CURSED_LOG.print(*args, **kwargs)



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



def cursed(fun, *args, **kwargs):

    def cursed_fun(*args, **kwargs):
        try:
            return fun(*args, **kwargs)
        except Exception as e:
            try:
                curses.nocbreak()
                curses.echo()
                curses.nonl()
                curses.endwin()
                log("warning", "Curse escaped")
            except Exception as ce:
                log("warning", "Curses did not close properly:", ce)
            log("warning", "Cursed error")
            raise e
    return cursed_fun


class CursedWindow(object):
    def __init__(self, parent, height=None, width=None, y=None, x=None, box=True, title=None):
        if height is None:
            height = curses.LINES
        if width is None:
            width = curses.COLS
        if x is None:
            x = 0
        if y is None:
            y = 0

        if type(width) is list:
            width = width[0]
        if type(height) is list:
            height = height[0]
        if type(x) is list:
            x = x[1]
        if type(y) is list:
            y = y[1]

        self.parent = parent
        self.height, self.width, self.y, self.x = height, width, y, x
        self.box = box

        self.window = curses.newwin(self.height, self.width, self.y, self.x)
        self.content = []

        self.colours = {
            "default": {"n": 0, "f": curses.COLOR_WHITE, "b":curses.COLOR_BLACK},
            "black": {"n": 1, "f": curses.COLOR_BLACK, "b":curses.COLOR_WHITE},
            "blue": {"n": 2, "f": curses.COLOR_BLUE, "b":curses.COLOR_WHITE},
            "cyan": {"n": 3, "f": curses.COLOR_CYAN, "b":curses.COLOR_BLACK},
            "green": {"n": 4, "f": curses.COLOR_GREEN, "b":curses.COLOR_BLACK},
            "magenta": {"n": 5, "f": curses.COLOR_MAGENTA, "b":curses.COLOR_BLACK},
            "red": {"n": 6, "f": curses.COLOR_RED, "b":curses.COLOR_BLACK},
            "white": {"n": 7, "f": curses.COLOR_WHITE, "b":curses.COLOR_BLACK},
            "yellow": {"n": 8, "f": curses.COLOR_YELLOW, "b":curses.COLOR_BLACK},
        }
        self.init_colours()


        if box:
            self.window.box()
            self.height -= 2
            self.width -= 2

        self.set_title(title)

        self.refresh()

    def refresh(self):
        self.window.refresh()

    @cursed
    def init_colours(self):
        for n, (k, v) in enumerate(self.colours.items()):
            curses.init_pair(n+1, v["f"], v["b"])
            v["c"] = n+1

    def get_colour(self, name):
        return self.colours.get(name, {}).get("c", 0)

    @cursed
    def print(self, *vals, c=None, end="\n"):
        val = " ".join(vals)
        for n, line in enumerate(val.split("\n")):
            if len(line) > self.width:
                line = line[:self.width]
            line = line

            line = {"l":line, "c":self.get_colour(c), "end":end}
            self.content.append(line)
        if len(self.content) > self.height:
            self.content = self.content[len(self.content)-self.height:]

        y = 1
        x = 1
        for n, s in enumerate(self.content):
            col = curses.color_pair(s.get("c", 0))
            
            extra_space = self.width-len(s["l"]) - (x-2)
            if extra_space < 0:
                l = s["l"][:extra_space]
            else:
                l = s["l"]+" ".join(["" for _ in range(extra_space)])
            try:
                self.window.addstr(y, x, l, col)
            except:
                raise Exception(x, y, l, col, extra_space)

            end = s["end"]
            if end == "\n":
                y += 1
                x = 1
            elif end == "\r":
                x = 1
            elif end == " ":
                x = x + len(s["l"])+1
            elif end == "":
                x = x + len(s["l"])
        self.refresh()

    @cursed
    def set_title(self, title=None):
        if title is not None:
            title = str(title)
            if len(title) > self.width-4:
                title = title[:self.width-4]
        if self.box:
            self.window.box()
            if title is not None:
                self.window.addstr(0, 2, f" {title} ")
        self.refresh()





class CursedTerminal(object):
    def __init__(self):
        
        self.screen = None
        self.windows = []


        self._init_screen()
    @cursed
    def _init_screen(self, screen=None):
            if screen is None:
                screen = curses.initscr()
                curses.start_color()
            self.screen = screen
    @cursed
    def close(self):
        curses.nocbreak()
        #self.screen.keypad(0)
        curses.echo()
        curses.nonl()
        for win in self.windows():
            win.endwin()
        curses.endwin()
        CURSED = False
        CURSED_LOG = None
        log("warning", "Curse closed")


    @cursed
    def refresh(self):
        self.screen.refresh()

    @cursed
    def refresh_windows(self):
        for window in self.windows.values():
            window.refresh()
    @cursed
    def add_window(self, *args, window_class=CursedWindow,**kwargs):
        w = window_class(self, *args, **kwargs)
        self.windows.append(w)
        return w
    @cursed
    def split_width(self, divs=2, split_height=False, percentages=None):
        if percentages is None:
            percentages = [100/divs]*divs
        #print(percentages)
        if split_height:
            w = curses.LINES
        else:
            w = curses.COLS 
        ws = [math.floor(w*p/100) for p in percentages]
        totalw = 0
        tw = []
        for ww in ws:
            tw.append([ww, totalw])
            totalw+=ww
        if split_height:
            tw[-1][0] += curses.LINES - totalw
        else:
            tw[-1][0] += curses.COLS - totalw
        return tw

    @cursed
    def split_height(self, *args, **kwargs):
        return self.split_width(*args, split_height=True, **kwargs)





class _deprecated():
    def update(self, val):
        from curses import wrapper
        import curses
        try:
            self.vals.append(val)
            max_vals = (os.get_terminal_size().columns // self.col_width) -2

            if len(self.vals) > max_vals:
                self.vals = self.vals[len(self.vals)-max_vals:]
            vals = np.array(self.vals)
            self.title.addstr(1, 1, str(vals))
            
            highest = max(vals)

            vals = (vals // (highest/(self.height)))
            #print(vals)
            #wrapper(self._update)
        
            self.title.addstr(2, 1, f"self.vals: {self.vals}")
            self.title.addstr(3, 1, f"vals:      {vals}")
            self.title.addstr(4, 1, f"highest: {highest}")
            self.title.refresh()
            self._update_graph(self.graph, vals)
        except:
            curses.nocbreak()
            self.screen.keypad(0)
            curses.echo()
            curses.endwin()
            curses.nonl()
            raise

    def _update_graph(self, graph, vals):

        #print(screen)
        for n, val in enumerate(vals):
            #print(n, val)
            for i in range(int(val)):
                #raise Exception(self.height-int(i)-2)
                graph.addstr(self.height-int(i), n+1, "#")
        graph.refresh()

def quad_term():
    term = CursedTerminal()

    heights = term.split_height(percentages=[10,50,40])
    #print(heights)
    widths = term.split_width(2)
    #print(widths)

    top =    term.add_window(heights[0], None,      heights[0], 0,         title=f"Top")
    left =   term.add_window(heights[1], widths[0], heights[1], widths[0], title=f"Left")
    right =  term.add_window(heights[1], widths[1], heights[1], widths[1], title=f"Right")
    bottom = term.add_window(heights[2], None,      heights[2], 0,         title=f"Bottom")

    return term, top, left, right, bottom