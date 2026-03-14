#! /usr/bin/env python
# -*- coding: utf-8 -*-

#future imports and compatibility between 2.x and 3.x
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

from future.standard_library import install_aliases
install_aliases()

import sys
info_p = sys.version_info
info_g = (sys.version).splitlines()
PYTHON_V = info_p.major

if PYTHON_V == 3:
    from builtins import input
    from builtins import object
    from builtins import range
    from builtins import str

from urllib.parse import urlparse, urlencode
from urllib.request import urlopen, Request, build_opener, install_opener, HTTPRedirectHandler, ProxyHandler
from urllib.error import HTTPError

import ctypes
import cProfile
import inspect
import math
import os
import platform
import signal
import subprocess
import threading
import time
import traceback
import re
import multiprocessing
import numpy
import matplotlib.pyplot as plt
import psutil
import io
import logging
import logging.handlers
logger = logging.getLogger("arcimboldo")

from termcolor import colored

#import managerSSH

DicGrid = None
cmLo = None
NODES = []
LISTJOBS = {}
list_commands = []


#######################################################################################################
#                                           SUPPORT FUNC                                              #
#######################################################################################################


string_types = (type(b''), type(u''))
import functools
import inspect
import warnings

def deprecated(reason):
    """
    This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used.
    """

    if isinstance(reason, string_types):

        # The @deprecated is used with a 'reason'.
        #
        # .. code-block:: python
        #
        #    @deprecated("please, use another function")
        #    def old_function(x, y):
        #      pass

        def decorator(func1):

            if inspect.isclass(func1):
                fmt1 = "Call to deprecated class {name} ({reason})."
            else:
                fmt1 = "Call to deprecated function {name} ({reason})."

            @functools.wraps(func1)
            def new_func1(*args, **kwargs):
                warnings.simplefilter('always', DeprecationWarning)
                warnings.warn(
                    fmt1.format(name=func1.__name__, reason=reason),
                    category=DeprecationWarning,
                    stacklevel=2
                )
                warnings.simplefilter('default', DeprecationWarning)
                return func1(*args, **kwargs)

            return new_func1

        return decorator

    elif inspect.isclass(reason) or inspect.isfunction(reason):

        # The @deprecated is used without any 'reason'.
        #
        # .. code-block:: python
        #
        #    @deprecated
        #    def old_function(x, y):
        #      pass

        func2 = reason

        if inspect.isclass(func2):
            fmt2 = "Call to deprecated class {name}."
        else:
            fmt2 = "Call to deprecated function {name}."

        @functools.wraps(func2)
        def new_func2(*args, **kwargs):
            warnings.simplefilter('always', DeprecationWarning)
            warnings.warn(
                fmt2.format(name=func2.__name__),
                category=DeprecationWarning,
                stacklevel=2
            )
            warnings.simplefilter('default', DeprecationWarning)
            return func2(*args, **kwargs)

        return new_func2

    else:
        raise TypeError(repr(type(reason)))

def timing(f):
    """

    :param f:
    :type f:
    :return:
    :rtype:
    """

    def wrap(*args, **kwds):
        time1 = time.time()
        ret = f(*args, **kwds)
        time2 = time.time()
        # if time2-time1>=1:
        print('%s function took %0.3f s' % (f.__name__, (time2 - time1)))
        return ret

    return wrap

def profileit(func):
    """

    :param func:
    :type func:
    :return:
    :rtype:
    """

    def wrapper(*args, **kwargs):
        datafn = func.__name__ + ".profile"  # Name the data file sensibly
        prof = cProfile.Profile()
        retval = prof.runcall(func, *args, **kwargs)
        prof.dump_stats(datafn)
        return retval

    return wrapper

def request_url(search_string, rootsite, data=None):
    z = None
    tries = 10
    for trial in range(tries):
        try:
            if data is None:
                #print("METHOD A", search_string, data)
                req = Request(search_string)
                z = urlopen(req).read()
            else:
                #print("METHOD B", search_string, data)
                req = Request(search_string)
                z = urlopen(req, data).read()
            break
        except:
            # print(sys.exc_info())
            # traceback.print_exc(file=sys.stdout)
            print("Error contacting",rootsite,"...Attempt ", trial + 1, "/", tries, search_string, data)
            #time.sleep(1)
    if not z:
        print("Connection Error")
        return None
    else:
        return z

def plot_scoring_fn(u, mean, topx, v=0.9, p=0.5, n=0.5, step=0.01, minx=-0.5, topxx=3.0, topy=1.2, miny=-0.2,
                    tickx=0.5, ticky=0.2):
    """

    :param u:
    :type u:
    :param mean:
    :type mean:
    :param topx:
    :type topx:
    :param v:
    :type v:
    :param p:
    :type p:
    :param n:
    :type n:
    :param step:
    :type step:
    :return:
    :rtype:
    """
    t = [(u, numpy.abs((u - mean) / mean) ** n if u <= mean else (((((u - mean) * ((mean * v))) / (
        mean + ((topx - mean) * p) - mean)) / mean)) ** n if u <= mean + ((topx - mean) * p) else (v + (((((u - (
        mean + (topx - mean) * p)) * (mean - (mean * v))) / (topx - (mean + (topx - mean) * p) + (
        mean * v))) / mean))) ** n) for u in numpy.arange(0, topx, step)]
    plt.xlim(minx, topxx)
    plt.ylim(miny, topy)
    plt.scatter(*list(zip(*t)), color='black')
    plt.show()
    print(t)
    # fig = plt.figure()
    # fig.savefig("./hola.png")
    return t

def safe_get_from_list(_list, index, default=None):
    try:
        return _list[index]
    except:
        return default

def safe_call(funct, *args, **kwargs):
    try:
        return funct(args, kwargs)
    except:
        return None

def py2_3_unicode(value):
    global PYTHON_V

    if PYTHON_V == 2:
        return unicode(value)
    else:
        return value

class FileNotFound(Exception): pass

class GiveMeNext(object):
    def __init__(self, iterable):
        self._iter = iter(iterable)

    def __next__(self):  # Py3-style iterator interface
        return next(self._iter)  # builtin next() function calls

    def __iter__(self):
        return self
#######################################################################################################
#                                                                                                     #
#######################################################################################################


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

def warning(objs):
    sys.stderr.write("WARNING: " + str(objs) + "\n")


def findInSubdirectory(filename, subdirectory=''):
    if subdirectory:
        path = subdirectory
    else:
        path = os.getcwd()
    for root, subFolders, files in os.walk(path):
        for fileu in files:
            if fileu == filename:
                return os.path.join(root, fileu)
    raise FileNotFound()


def startCheckQueue(sym, delete_check_file=False, callback=None, forcework=False):
    """ This function allows to create a queue when you do not have a queue manager.

    It is used in the ARCIMBOLDO_* programs in the supercomputer mode, and in BORGES_MATRIX to parallelize functions

    :param sym:
    :type sym:
    :param delete_check_file:
    :type delete_check_file:
    :param callback:
    :type callback:
    :param forcework:
    :type forcework:
    :return:
    :rtype:
    """
    global NODES

    if not forcework and len(NODES) == 0:
        return

    if sym.PROCESSES > 0:
        print("I found ", sym.REALPROCESSES, "CPUs.")
        while 1:
            if len(threading.enumerate()) <= sym.PROCESSES:
                p = OutputThreading(__startQueue, delete_check_file, callback)
                p.start()
                break
    else:
        print("FATAL ERROR: I cannot load correctly information of CPUs.")
        # sym.couldIClose = True
        sys.exit(0)


def checkYOURoutput(myfile, conditioEND, testEND, sleep_ifnot_ready=True, failure_test=None):
    esegui = True
    correct = False
    while esegui:
        if not os.path.exists(myfile):
            if sleep_ifnot_ready:
                time.sleep(3)
                continue
            else:
                return False
        if failure_test != None and testEND != None and sleep_ifnot_ready:
            f = open(myfile)
            e = f.read()
            f.close()
            if isinstance(failure_test, str):
                numero = e.count(failure_test)
                # print "--------",myfile,failure_test,numero,testEND,conditioEND
                if numero == int(testEND):
                    return None
            elif isinstance(conditioEND, list):
                numero = 0
                for uu in conditioEND:
                    numero += e.count(uu)
                if numero == int(testEND):
                    return None
        if conditioEND != None and testEND != None:
            f = open(myfile)
            e = f.read()
            f.close()
            if isinstance(conditioEND, str):
                numero = e.count(conditioEND)
                # print "+++--------+++",myfile,numero,conditioEND,testEND,type(numero),type(testEND)
                if numero == int(testEND):
                    esegui = False
                    correct = True
                elif sleep_ifnot_ready:
                    time.sleep(3)
                else:
                    return False
            elif isinstance(conditioEND, list):
                numero = 0
                for uu in conditioEND:
                    numero += e.count(uu)
                if numero == int(testEND):
                    esegui = False
                    correct = True
                elif sleep_ifnot_ready:
                    time.sleep(3)
                else:
                    return False
        else:
            esegui = False
            correct = True
    return correct


def __startQueue(delete_check_file=False, callback=None):
    global NODES
    global LISTJOBS

    while 1:
            node = None
            if len(LISTJOBS.keys()) == 0: continue
            else: node = list(LISTJOBS.keys())[0]
            joblist = LISTJOBS[node]
            test = False
            done = []
            for j in range(len(joblist)):
                job, conditio, end_number, blocking, command = joblist[j]
                # print "Checking job",job,conditio,end_number
                if os.path.exists(job):
                    test = checkYOURoutput(job, conditio, end_number, sleep_ifnot_ready=False)
                    # print "Testing ",job,"test",test,conditio,end_number
                if test:
                    done.append(j)
                    if delete_check_file:
                        os.remove(job)
                    if callback != None:
                        callback(os.path.basename(job)[:-4])
            # print "len(done) is",len(done),"len(joblist) is",len(joblist)
            if len(done) == len(joblist):
                print("Freeing the node:", node, "from jobs:", len(done))
                del LISTJOBS[node]
            else:
                LISTJOBS[node] = multi_delete(LISTJOBS[node], done)

                # for j in done:
                #       del LISTJOBS[node][j]

def endCheckQueue(blocking=True):
    global list_commands
    global LISTJOBS

    if len(list_commands) > 0:
        valus = list_commands[-1]
        del list_commands[-1]
        launchCommand(valus[0], valus[1], valus[2], valus[3], single=False, sendnow=True, blocking=valus[4])

    if blocking:
        while len(LISTJOBS.keys()) > 0:
            # print "Waiting for",len(LISTJOBS.keys()),"to end..."
            time.sleep(3)

        print ("Queue is Empty!")
    else:
        print ("Not blocking queue check")


def slice_list(inputl, size):
    input_size = len(inputl)
    slice_size = input_size / size
    remain = input_size % size
    result = []
    iterator = iter(inputl)
    for i in range(size):
        result.append([])
        for j in range(slice_size):
            result[i].append(next(iterator))
        if remain:
            result[i].append(next(iterator))
            remain -= 1
    return result


def __sendJob(dicjob):
    for nodo in sorted(dicjob.keys()):
        comandi = dicjob[nodo]
        for comm in comandi:
            times = 1
            while 1:
                try:
                    conn = managerSSH.Connection(nodo)
                    break
                except:
                    # print sys.exc_info()
                    # traceback.print_exc(file=sys.stdout)
                    if times <= 10:
                        times += 1
                        time.sleep(1)
                    else:
                        break
            # channel = conn.interactive()
            # print "nodo",nodo,"comando start...",conn
            # print "",comm
            # conn.send_command_to_channel(channel,comm,"$ ")
            print (conn.execute(comm))
            # sys.stdout.readlines()
            print ("comando end...")
            # channel.close()
            conn.close()

def launchCommand(command, valore, conditio, end_number, single=False, sendnow=False, blocking=True):
    global NODES
    global LISTJOBS
    global list_commands

    list_commands.append((command, valore, conditio, end_number, blocking))
    if not sendnow and not single:
        return

    # print "NODES are:========="
    # print NODES
    # print "==================="
    emptyNODES = []
    # eventlet.monkey_patch(os=False, thread=False)
    while 1:
        if len(LISTJOBS.keys()) < len(NODES):
            found = False
            nodo = ""
            ns = ""
            for n in NODES:
                if n not in LISTJOBS.keys():
                    found = True
                    nodo = n.split("***")[0]
                    ns = n
                    # print "Nodo selected is",nodo,ns,"Nodes occupied",len(LISTJOBS.keys())
                    emptyNODES.append((nodo, n))

            if not found:
                print ("All nodes are occupied", len(LISTJOBS))
                time.sleep(1)
                continue

            print ("Executing...")
            # user = getpass.getuser()
            # first = "ssh "+user+"@"+nodo+' '+'"cd '+os.path.split(valore)[0]+';nohup '+command+' &" &'
            if not single:
                listas = slice_list(list_commands, len(emptyNODES))
            else:
                listas = [list_commands]

            hwmany = {}
            for lis in range(len(listas)):
                list_valo = []
                first = ""
                lisi = listas[lis]
                nodo = emptyNODES[lis][0]
                n = emptyNODES[lis][1]
                # print "Sending",len(lisi),"jobs in the core",n
                for it in range(len(lisi)):
                    command, valore, conditio, end_number, blocking = lisi[it]
                    list_valo.append((valore, conditio, end_number, blocking, command))
                    if os.path.exists(os.path.split(valore)[0]):
                        first += 'cd ' + os.path.split(valore)[
                            0] + ';nohup ' + command  # +' 2> '+os.path.join(os.path.split(valore)[0],'./error.txt') #/dev/null'
                    else:
                        first += 'nohup ' + command  # +' 2> /dev/null'
                    if it < len(lisi) - 1:
                        first += " ; "
                    else:
                        first += " & "
                # print first
                if nodo in hwmany.keys():
                    hwmany[nodo].append(first)
                else:
                    hwmany[nodo] = [first]

                LISTJOBS[n] = list_valo
            __sendJob(hwmany)
            list_commands = []
            break
        else:
            # print "All nodes are occupied, sorry",LISTJOBS
            time.sleep(2)
            # eventlet.monkey_patch(all=True)


def floatRgb(mag, cmin, cmax):
    """
    Return a tuple of floats between 0 and 1 for the red, green and
    blue amplitudes.
    """

    try:
        # normalize to [0,1]
        x = float(mag - cmin) / float(cmax - cmin)
    except:
        # cmax = cmin
        x = 0.5
    blue = min((max((4 * (0.75 - x), 0.)), 1.))
    red = min((max((4 * (x - 0.25), 0.)), 1.))
    green = min((max((4 * math.fabs(x - 0.5) - 1., 0.)), 1.))
    return (red, green, blue)


def strRgb(mag, cmin, cmax):
    """
    Return a tuple of strings to be used in Tk plots.
    """

    red, green, blue = floatRgb(mag, cmin, cmax)
    return "#%02x%02x%02x" % (red * 255, green * 255, blue * 255)


def rgb(mag, cmin, cmax):
    """
    Return a tuple of integers to be used in AWT/Java plots.
    """

    red, green, blue = floatRgb(mag, cmin, cmax)
    return (int(red * 255), int(green * 255), int(blue * 255))


def htmlRgb(mag, cmin, cmax, tone):
    """
    Return a tuple of strings to be used in HTML documents.
    """
    try:
        # normalize to [0,1]
        x = float(mag - cmin) / float(cmax - cmin)
    except:
        # cmax = cmin
        x = 0.5

    if tone == "red":
        red = int(255 * x)
        green = int(255 * (1.0 - x))
        blue = 0
    elif tone == "green":
        red = 0
        green = int(255 * x)
        blue = int(255 * (1.0 - x))
    else:
        red = int(255 * (1.0 - x))
        green = 0
        blue = int(255 * x)

    return "#%02x%02x%02x" % (red, green, blue)


# TODO: Comment all the functions of this module. Claudia will search the correct way of doing it.
def open_connection(DicGridConn, DicParameters, cm):
    """
    open_connection(dict:DicGridConn, dict:DicParameters, Grid:cm)
    **Open the remote connection with an host
    **Store DicGridConn and cm as global parameter
    **If DicGridConn is empty does nothing
    ============================
    DicGridConn:
            ---key : value---
            username: string. the username
            host: string. Is the hostname
            port: integer. Is the port to connect to the remote host
            passkey: string. Could be a password or a passkey file path or the string False
            promptA: prompt linux shell of frontend
            remote_submitter_password: string. Password for the submitter host
            remote_submitter_username:
            remote_submitter_host:
            remote_submitter_port:
            home_frontend_directory:
            promptB: prompt linux shell of submitter
            isnfs: boolean. It specifies wether the filesystem into the remote host in nfs or not.
    DicParameters:
            key : value
    cm:
            Grid Object that manages the middleware grid into the remote host
    =============================
    Return:
         void
    """

    global DicGrid
    global cmLo

    # TODO: Transform this module or part of this module in a class to avoid global parameters and convert it
    #      int attribute class. Draw the Class.
    if len(DicGridConn.keys()) > 0:
        DicGridConn, cm = requestConnectionRemote(DicGridConn, cm)
        DicGrid = DicGridConn
        cmLo = cm
        print (cm.create_remote_dir(DicParameters["nameExecution"]))
        print ("During opening I tried to create", DicParameters["nameExecution"])
        print (cm.change_remote_dir(DicParameters["nameExecution"]))


def close_connection(DicGridConn, DicParameters, cm):
    if hasattr(cm, "channel"):
        print (cm.change_remote_dir(".."))
        print (cm.remove_remote_dir(DicParameters["nameExecution"]))
        print (cm.close_exchange_file())
        print (cm.close_remote_connection())


def remote_reconnection(directory):
    global DicGrid
    global cmLo

    cmLo.lock.acquire()
    print (cmLo.close_exchange_file())
    print (cmLo.close_remote_connection())

    if len(DicGrid.keys()) > 0:
        DicGrid, cmLo = requestConnectionRemote(DicGrid, cmLo)
        print (cmLo.change_remote_dir(directory))
    cmLo.lock.release()


def requestConnectionRemote(DicGridConn, cm):
    times = 10
    count = 0

    while True:
        try:
            command = ""
            if DicGridConn["passkey"] == "False":
                # TODO: Keywords of DicGridConn should be consistent in names
                print ("Please insert the password for ", DicGridConn["username"] + "@" + DicGridConn[
                    "host"], "with port", DicGridConn["port"])
                command = input("<> ")
                assert isinstance(command, str)

                # NOTE: Do we really want reuse the same key for a different meaning?
                DicGridConn["passkey"] = {"password": command}
                command = ""
            listona = []
            if "remote_submitter_username" in DicGridConn.keys() and "remote_submitter_host" in DicGridConn.keys() and "promptB" in DicGridConn.keys() and \
                                            DicGridConn["remote_submitter_username"] + "@" + DicGridConn[
                                "remote_submitter_host"] != DicGridConn["username"] + "@" + DicGridConn["host"]:
                if "remote_submitter_password" not in DicGridConn.keys():
                    if "remote_submitter_port" not in DicGridConn.keys():
                        DicGridConn["remote_submitter_port"] = 22
                    print ("Please insert the password for ", DicGridConn["remote_submitter_username"] + "@" + \
                                                             DicGridConn["remote_submitter_host"], "with port", \
                    DicGridConn["remote_submitter_port"])
                    command = input("<> ")
                    assert isinstance(command, str)
                else:
                    command = DicGridConn["remote_submitter_password"]

                DicGridConn["remote_submitter_password"] = command
                listona = [
                    ["sublocalssh", DicGridConn["remote_submitter_username"], DicGridConn["remote_submitter_host"],
                     DicGridConn["remote_submitter_port"], DicGridConn["remote_submitter_password"],
                     DicGridConn["promptB"]]]

            amm = DicGridConn["home_frontend_directory"]
            if amm.strip() in ["", " ", None]:
                amm = "./"

            if "promptB" not in DicGridConn.keys():
                DicGridConn["promptB"] = DicGridConn["promptA"]

            DicGridConn["listInitialOp"] = listona
            out = cm.getGridJobsToRemoteMachine(DicGridConn["username"], DicGridConn["host"], DicGridConn["port"],
                                                DicGridConn["passkey"], DicGridConn["promptA"], DicGridConn["promptB"],
                                                DicGridConn["isnfs"],
                                                listOfInitialOperations=DicGridConn["listInitialOp"],
                                                home_frontend_directory=amm)
            if not out.strip().endswith('password:'):
                return DicGridConn, cm
            else:
                cm.close_exchange_file()
                cm.close_remote_connection()
                #cm = None
                #cm = Grid.condorManager()
                raise Exception()
        except:
            count += 1
            print ("An error occured while connecting on remote machine PID: " + str(os.getpid()) + " " + str(
                count) + " time.")
            if count >= times:
                if DicGridConn["passkey"] != "False":
                    DicGridConn["passkey"] = "False"
                    count = 0
                    continue
                else:
                    print (sys.exc_info())
                    traceback.print_exc(file=sys.stdout)
                    # print paramiko.BadAuthenticationType.allowed_types
                    couldIClose = True
                    time.sleep(2)
                    sys.exit(1)


# def multi_delete(list_, *args):
#     indexes = sorted(list(args), reverse=True)
#     for index in indexes:
#         del list_[index]
#     return list_

def multi_delete(list_, args):
    indexes = sorted(args, reverse=True)
    for index in indexes:
        del list_[index]
    return list_


def mem(size="rss"):
    """Generalization; memory sizes: rss, rsz, vsz."""
    return int(os.popen('ps -p %d -o %s | tail -1' % (os.getpid(), size)).read())


def rss():
    """Return ps -o rss (resident) memory in kB."""
    return mem("rss")


def rsz():
    """Return ps -o rsz (resident + text) memory in kB."""
    return mem("rsz")


def vsz():
    """Return ps -o vsz (virtual) memory in kB."""
    return mem("vsz")

def generateLogger():
    """
    Create logger: In case there is no working directory, the information
    will be stored in a buffer instead of a file. The buffer can be dumped to
    a file later.
    """
    logstream = io.StringIO()

    streamHandler = logging.StreamHandler(logstream)
    streamHandler.setLevel(logging.DEBUG)
    logger.addHandler(streamHandler)

    stdoutHandler = logging.StreamHandler(sys.stdout)
    stdoutHandler.setLevel(logging.DEBUG)
    logger.addHandler(stdoutHandler)

    logger.setLevel(logging.DEBUG)

def generateLoggerDirectory(workingDirectory, logName="arcimboldo.log"):
    """
    Create logger in a working directory with a specific name:
        - Check if this function has already been called, if so, return
        - Remove previous logger handlers
        - Rotate files in the working directory
        - Maximum size of 1GB
        - logger.debug: Store information in the logging file
        - logger.info: Stdout and store information in the logging file
        - Move previous information from generateLogger to the file
    """

    loggerExist = False
    loggerData = ""
    for handler in logger.handlers:  # check if exists the logger already
        if isinstance(handler, logging.FileHandler):
            return
        loggerExist = True

    if loggerExist:
        loggerData = logger.handlers[0].stream.getvalue()
        for handler in logger.handlers[:]:  # remove all old handlers
            logger.removeHandler(handler)

    logPath = os.path.join(workingDirectory, logName)
    handler = logging.handlers.RotatingFileHandler(logPath, mode='w', maxBytes=1 * 1024 * 1024 * 1024, backupCount=5, delay=True)
    if os.path.isfile(logPath):
        handler.doRollover()
    logger.setLevel(logging.DEBUG)

    infoHandler = logging.StreamHandler(sys.stdout)
    infoHandler.setLevel(logging.INFO)
    logger.addHandler(infoHandler)

    debugHandler = logging.FileHandler(logPath)
    debugHandler.setLevel(logging.DEBUG)
    logger.addHandler(debugHandler)

    if loggerData != "":
        logger.debug(loggerData)

# TODO: Managing Threads should be always in synchro. the main thread should have the possibility to know
# when a thread child has died, finished or exit and take decisions by consequence.

class OutputThreading(threading.Thread):
    def __init__(self, externalCallable, *args, **kwds):
        threading.Thread.__init__(self, **kwds)
        self.externalCallable = externalCallable
        self.args = args
        self.kwds = kwds

    def run(self):
        ar = self.args
        kw = self.kwds
        self.externalCallable(*ar, **kw)

    def _get_my_tid(self):
        """determines this (self's) thread id

        CAREFUL : this function is executed in the context of the caller
        thread, to get the identity of the thread represented by this
        instance.
        """
        if not self.isAlive():
            raise threading.ThreadError("the thread is not active")

        # do we have it cached?
        if hasattr(self, "_thread_id"):
            return self._thread_id

        # no, look for it in the _active dict
        for tid, tobj in threading._active.items():
            if tobj is self:
                self._thread_id = tid
                return tid

        # TODO: in python 2.6, there's a simpler way to do : self.ident

        raise AssertionError("could not determine the thread's id")

    def raiseExc(self, exctype):
        """Raises the given exception type in the context of this thread.

        If the thread is busy in a system call (time.sleep(),
        socket.accept(), ...), the exception is simply ignored.

        If you are sure that your exception should terminate the thread,
        one way to ensure that it works is:

            t = ThreadWithExc( ... )
            ...
            t.raiseExc( SomeException )
            while t.isAlive():
                time.sleep( 0.1 )
                t.raiseExc( SomeException )

        If the exception is to be caught by the thread, you need a way to
        check that your thread has caught it.

        CAREFUL : this function is executed in the context of the
        caller thread, to raise an excpetion in the context of the
        thread represented by this instance.
        """
        _async_raise(self._get_my_tid(), exctype)


def _async_raise(tid, exctype):
    '''Raises an exception in the threads with id tid'''
    if not inspect.isclass(exctype):
        raise TypeError("Only types can be raised (not instances)")
    res = ctypes.pythonapi.PyThreadState_SetAsyncExc(ctypes.c_long(tid),
                                                     ctypes.py_object(exctype))
    if res == 0:
        raise ValueError("invalid thread id")
    elif res != 1:
        # "if it returns a number greater than one, you're in trouble,
        # and you should call it again with exc=NULL to revert the effect"
        ctypes.pythonapi.PyThreadState_SetAsyncExc(ctypes.c_long(tid), None)
        raise SystemError("PyThreadState_SetAsyncExc failed")

def signal_term_handler(signal,frame):
    print ('\nThe program received a signal of termination')
    if hasattr(sys, '_MEIPASS'):
        print ("Temp file was ", sys._MEIPASS, " and was removed before exiting")
    sys.exit(0)


class SystemUtility:
    PROCESSES = 0
    couldIClose = False
    REALPROCESSES = 0

    def __init__(self,hyper=False):
        self.hyper = hyper
        if not self.hyper:
            number_cpus = psutil.cpu_count(logical = False) - 1
            if number_cpus is None:
                number_cpus = self.cpu_count() - 1
            self.PROCESSES = number_cpus
            self.REALPROCESSES = self.PROCESSES
        else:
            self.PROCESSES = multiprocessing.cpu_count() - 1
            self.REALPROCESSES = self.PROCESSES
        signal.signal(signal.SIGUSR1, self.incrementProcesses)
        signal.signal(signal.SIGUSR2, self.decrementProcesses)

    def startCheckingSystem(self):
        OutputThreading(self.checkMemoryAndLAUsed).start()
        pass

    def incrementProcesses(self, signal, frame):
        print ('You ask to increment the number of processes!')
        self.PROCESSES += 1
        if self.PROCESSES <= 3:
            self.PROCESSES = 3
        elif self.PROCESSES >= 200:
            self.PROCESSES = 200

    def decrementProcesses(self, signal, frame):
        print ('You ask to decrement the number of processes!')
        self.PROCESSES -= 1
        if self.PROCESSES <= 3:
            self.PROCESSES = 3
        elif self.PROCESSES >= 200:
            self.PROCESSES = 200

    def checkMemoryAndLAUsed(self):
        while not self.couldIClose:
            oldValue = self.PROCESSES
            # NOTE:LINUX
            if platform.system().lower() == "linux":
                # p0 = subprocess.Popen(["flush"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                p1 = subprocess.Popen(["free"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                p2 = subprocess.Popen(["grep", "Mem"], stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
                out, err = p2.communicate()
                liMem = out.split()
                total = int(liMem[1])
                used = int(liMem[2])
                free = int(liMem[3])
                buffers = int(liMem[5])
                cached = int(liMem[6])
            elif platform.system().lower() in ["windows", "win32"]:
                # NOTE: WINDOZZZ
                import psutil
                tupn = psutil.virtual_memory()
                total = tupn.total / 1024.0
                free = tupn.available / 1024.0
            elif platform.system().lower() == "darwin":
                # NOTE:MAC
                p1 = subprocess.Popen(["vm_stat"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                vm, err = p1.communicate()
                p1 = subprocess.Popen(["sysctl", "-n", "hw.memsize"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                outspl, err = p1.communicate()
                outspl = out.splitlines()
                for line in outspl:
                    line = line.strip()
                    try:
                        total = int(line.split()[0]) / 1024.0
                    except:
                        continue

                # Process vm_stat
                vmLines = vm.split('\n')
                sep = re.compile(':[\s]+')
                for row in range(1, len(vmLines) - 2):
                    rowText = vmLines[row].strip()
                    rowElements = sep.split(rowText)
                    if "Pages free" == rowElements[0]:
                        free = int(rowElements[1].strip('\.')) * 4096
            else:
                print ("FATAL ERROR: Your operating system is: ", platform.system(), platform.release(), ", which does not seem to be a standard Linux, Mac OSX or Windows")
                sys.exit(1)

            percent = -1
            if total >= 5000000:
                percent = 10
            else:
                percent = 20

            kilobytes = (total * percent) / 100
            # print "free:",free,"compared",kilobytes,"oldValue",oldValue, (free-kilobytes),free < kilobytes
            # NOTE: LINUX
            condition = None
            if platform.system().lower() == "linux":
                condition = (free + buffers + cached) < kilobytes
            # NOTE: WINDOZZZ
            elif platform.system().lower() in ["windows", "win32"]:
                condition = free < kilobytes
            # NOTE: MAC
            elif platform.system().lower() == "darwin":
                condition = free < kilobytes
            else:
                print ("FATAL ERROR: Your operating system is: ", platform.system(), platform.release(), ", which does not seem to be a standard Linux or Mac OSX")
                sys.exit(1)

            if condition:
                #print "low memory free: i'm reducing the number of cuncurrent processes"
                self.REALPROCESSES = 1
            else:
                self.REALPROCESSES = oldValue
            # print "valore di real processes",self.REALPROCESSES

            # load average automatic guardian
            # NOTE:LINUX
            if platform.system().lower() == "linux":
                p = subprocess.Popen(["top", "-n", "1"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # NOTE:WINDOWZZZ
            elif platform.system().lower() in ["windows", "win32"]:
                pass
            elif platform.system().lower() == "darwin":
                # NOTE:MAC
                p = subprocess.Popen(["top", "-n1", "-l1"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            else:
                print ("FATAL ERROR: Your operating system is: ", platform.system(), platform.release(), ", which does not seem to be a standard Linux or Mac OSX")
                sys.exit(1)

            out, err = p.communicate()
            lilv = out.split()
            load_average = -1
            for u in range(len(lilv)):
                # NOTE:LINUX
                condition2 = None
                if platform.system().lower() == "linux":
                    condition2 = lilv[u] == "load" and lilv[u + 1] == "average:"
                # NOTE:MAC
                elif platform.system().lower() not in ["windows", "win32"]:
                    condition2 = lilv[u] == "Load" and lilv[u + 1] == "Avg:"
                else:
                    print ("FATAL ERROR: Your operating system is: ", platform.system(), platform.release(), ", which does not seem to be a standard Linux or Mac OSX")
                    sys.exit(1)

                if condition2:
                    load_average = float((lilv[u + 2])[:-1])
            if load_average != -1 and load_average <= (self.PROCESSES + 0.5):
                self.REALPROCESSES = 1
            else:
                self.REALPROCESSES = oldValue

            time.sleep(5)
        print ("The Borges Guardian Deamon is now stopped.")

    def spawn_function_with_multiprocessing(self,target, args):
        try:
            if self.PROCESSES > 0:
                print("I found ", self.REALPROCESSES, "CPUs.")
                while 1:

                    if len(multiprocessing.active_children()) < self.PROCESSES:
                        p2 = multiprocessing.Process(target=target, args=args)
                        p2.start()
                        return p2
            else:
                print("FATAL ERROR: I cannot load correctly information of CPUs.")
                self.couldIClose = True
                sys.exit(0)
        except KeyboardInterrupt:
            print("The user requires to exit from the program.")
            self.couldIClose = True
            sys.exit(0)

    def cpu_count(self):
        cpuDict = {}
        # NOTE:LINUX
        if platform.system().lower() == "linux":
            try:
                #p = subprocess.Popen(["grep", "-A", "3", "physical", "/proc/cpuinfo"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                #out, err = p.communicate()
                outlines = ""
                out = ""
                with open("/proc/cpuinfo","r") as rdo:
                    outlines = rdo.readlines()
                for zd,linel in enumerate(outlines):
                    if "physical id" in linel:
                        out += linel
                        out += outlines[zd+1]
                        out += outlines[zd+2]
                        out += outlines[zd+3]
                out += "--\n"

                licpu = out.split("--")
                for item in licpu:
                    liCampi = item.split("\n")
                    physical_id = -1
                    cpu_cores = -1
                    for it2 in range(len(liCampi)):
                        item2 = liCampi[it2]
                        if item2 == "":
                            continue
                        entry = item2.split()
                        if entry[0] == "physical" and entry[1] == "id":
                            physical_id = int(entry[3])
                        if entry[0] == "cpu" and entry[1] == "cores":
                            cpu_cores = int(entry[3])
                    if physical_id != -1 and cpu_cores != -1 and physical_id not in cpuDict.keys():
                        cpuDict[physical_id] = cpu_cores
                totalcpu = 0
                for key in cpuDict.keys():
                    totalcpu += cpuDict[key]

                if totalcpu == 0:
                    raise Exception()
            except:
                #print sys.exc_info()
                traceback.print_exc(file=sys.stdout)
                p = subprocess.Popen(["grep", "-c", "processor", "/proc/cpuinfo"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = p.communicate()
                totalcpu = int(out.strip())

        # NOTE:MAC
        elif platform.system().lower() not in ["windows", "win32"]:
            p = subprocess.Popen(["sysctl", "-n", "hw.ncpu"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            outspl = out.splitlines()
            for line in outspl:
                line = line.strip()
                try:
                    totalcpu = int(line.split()[0])
                except:
                    continue
        else:
            print ("FATAL ERROR: Your operating system is: ", platform.system(), platform.release(), ", which does not seem to be a standard Linux or Mac OSX")
            sys.exit(1)

        return totalcpu


class ThreadSafeObject:
    """
    A class that makes any object thread safe.
    """

    def __init__(self, obj):
        """
        Initialize the class with the object to make thread safe.
        """
        self.lock = threading.RLock()
        self.object = obj

    def __getattr__(self, attr):
        self.lock.acquire()

        def _proxy(*args, **kargs):
            self.lock.acquire()
            answer = getattr(self.object, attr)(*args, **kargs)
            self.lock.release()
            return answer

        return _proxy


class Singleton:
    """
    A non-thread-safe helper class to ease implementing singletons.
    This should be used as a decorator -- not a metaclass -- to the
    class that should be a singleton.

    The decorated class can define one `__init__` function that
    takes only the `self` argument. Other than that, there are
    no restrictions that apply to the decorated class.

    To get the singleton instance, use the `Instance` method. Trying
    to use `__call__` will result in a `TypeError` being raised.

    The decorated class cannot be inherited from.

    """

    def __init__(self, decorated):
        self._decorated = decorated

    def Instance(self, *args, **kwargs):
        """
        Returns the singleton instance. Upon its first call, it creates a
        new instance of the decorated class and calls its `__init__` method.
        On all subsequent calls, the already created instance is returned.

        """
        try:
            return self._instance
        except AttributeError:
            self._instance = self._decorated(*args, **kwargs)
            return self._instance

    def __call__(self):
        """
        Call method that raises an exception in order to prevent creation
        of multiple instances of the singleton. The `Instance` method should
        be used instead.

        """
        raise TypeError('Singletons must be accessed through the `Instance` method.')

class File(object):
    def __init__(self, path):
        self.name = path

# TODO:  Test main to test and check each function

