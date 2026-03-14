#! /usr/bin/env python

#! /usr/bin/env python
# -*- coding: utf-8 -*-

#future imports and compatibility between 2.x and 3.x
from __future__ import print_function
from __future__ import division
#from __future__ import unicode_literals

import sys
info_p = sys.version_info
info_g = (sys.version).splitlines()
PYTHON_V =  info_p.major

if PYTHON_V == 3:
    from builtins import input
    from builtins import range
    from builtins import str

import copy
import datetime
import getpass
import os
import subprocess
import tarfile
import threading
import time as timep
import traceback
import re


from . import SystemUtility
#import managerSSH

# NOTE: Please remember that you really need to do not use spaces with the redirect commands > and < for SGE and MOAB otherwise managerSSH could confuse redirects
# with shell prompt

PATH_LOCAL_PYTHON_INTERPRETER = "/usr/bin/python"
PATH_LOCAL_PYTHON_INTERPRETER = "/usr/bin/python"
PATH_REMOTE_PYTHON_INTERPRETER = "/usr/bin/python"
PATH_REMOTE_SGEPY = ""


class Grid(object):
    """SuperClass of all the Grid Managers """
    JOBS = {}  # jobid:{timeStart:time,cmd:string,job:Job,gridCluster:nCluster,queue:nQueue}
    JOBS_DONE = {}
    remote_library_path = ""
    remote_mtz_path = ""
    remote_hkl_path = ""
    remote_mtzP1_path = ""
    remote_tncs_path = ""
    remote_norm_path = ""
    FILE_TO_COPY = []
    JOB_TO_SUBMIT = []
    CUMULATIVE_TRANSFERING = 1000
    EXCEPTION_TO_CUMULATIVE = [".mtz", ".hkl", ".tncs", ".norm", ".cmd", "r.gz"]
    FILE_RECEIVED = {}
    actualDirectory = ""
    type_grid = "abstract"

    def __init__(self):
        self.isExchangeFileOpen = False
        self.exchangeConn = None
        self.promptSFTP = "> "
        self.lo = threading.RLock()
        # self.lock = ADT.LogLock("Grid.py",threading.Condition(self.lo))
        self.lock = threading.Condition(self.lo)

    def getGridJobsToRemoteMachine(self, username, host, port, passkey, prompta, promptb, isnfs,
                                   listOfInitialOperations=[], home_frontend_directory="./"):
        if isinstance(passkey, dict):
            self.connection = managerSSH.Connection(host, username=username, password=passkey["password"], port=port)
        else:
            print("Trying connecting with passkey", passkey)
            self.connection = managerSSH.Connection(host, username=username, private_key=passkey, port=port)

        self.promptA = prompta
        self.promptB = promptb
        self.username = username
        self.host = host
        self.port = port
        self.passkey = passkey
        self.listInitial = listOfInitialOperations
        self.isnfs = isnfs
        self.home_directory = home_frontend_directory

        buff = ""
        if len(listOfInitialOperations) == 0:
            self.channel = self.connection.interactive()
            # return ""
        else:
            for opera in listOfInitialOperations:
                if "sublocalssh" == opera[0]:
                    self.channel, buffo = self.connection.open_sublocal_interactive_ssh(opera[1], opera[2], opera[3],
                                                                                        opera[4], opera[5])
                    buff += buffo
                else:
                    if hasattr(self, "channel"):
                        buff += self.connection.send_command_to_channel(self.channel, opera[0], opera[1])
                    else:
                        self.channel = self.connection.interactive()
                        buff += self.connection.send_command_to_channel(self.channel, opera[0], opera[1])
        buff += "\n"
        buff += self.open_exchange_file()
        return buff

    def open_exchange_file(self):
        buff = ""
        if not self.isExchangeFileOpen and self.exchangeConn == None:
            try:
                if self.isnfs:
                    self.exchangeConn = self.connection.open_sftp_channel()
                    self.isExchangeFileOpen = True
                    for opera in self.listInitial:
                        comma = opera[0].split()
                        if len(comma) == 2 and comma[0] == "cd":
                            self.change_remote_dir(comma[1])
                else:
                    self.exchangeConn = self.connection.open_sftp_channel()
                    self.sftpsubcha = self.connection.interactive()
                    buff += self.connection.send_command_to_channel(self.sftpsubcha,
                                                                    "rm -rf " + os.path.join(self.home_directory,
                                                                                             "borges_middle"),
                                                                    self.promptA)

                    self.exchangeConn.mkdir(os.path.join(self.home_directory, "borges_middle"))
                    self.exchangeConn.chdir(os.path.join(self.home_directory, "borges_middle"))
                    buff += self.connection.send_command_to_channel(self.sftpsubcha,
                                                                    "cd " + os.path.join(self.home_directory,
                                                                                         "borges_middle"), self.promptA)
                    for opera in self.listInitial:
                        comma = opera[0].split()
                        if "sublocalssh" == opera[0]:
                            bluff = self.connection.send_command_to_channel(self.sftpsubcha,
                                                                            "sftp -oPort=" + str(opera[3]) + " " + str(
                                                                                opera[1]) + "@" + str(opera[2]),
                                                                            '\'s password: ', stopif=self.promptSFTP)
                            if not isinstance(bluff, tuple):
                                buff += bluff
                                buff += self.connection.send_command_to_channel(self.sftpsubcha, opera[4],
                                                                                self.promptSFTP)
                            else:
                                buff += bluff[0]
                        elif len(comma) == 2 and comma[0] == "cd":
                            self.change_remote_dir(comma[1])

                        self.isExchangeFileOpen = True

                return buff
            except:
                print ("Error on opening sftp connection")
                # print sys.exc_info()
                # traceback.print_exc()
        else:
            return "File Exchanger Protocol is already open!"

    def close_remote_connection(self):
        try:
            self.connection.close()
            print ("Remote Connection Protocol closed!")
        except:
            print ("Error on closing ssh connection")
            # print sys.exc_info()
            # traceback.print_exc(file=sys.stdout)

    def close_exchange_file(self):
        buff = ""
        if self.isExchangeFileOpen and self.exchangeConn != None:
            try:
                if self.isnfs:
                    self.exchangeConn.close()
                    self.isExchangeFileOpen = False
                else:
                    buff += self.connection.send_command_to_channel(self.sftpsubcha, "exit", self.promptA)
                    self.sftpsubcha.close()
                    self.exchangeConn.chdir("..")
                    # print self.exchangeConn.getcwd()
                    self.exchangeConn.rmdir("./borges_middle")
                    self.exchangeConn.close()
                    self.isExchangeFileOpen = False
                print ("File Exchanger Protocol closed!\n" + buff)
                self.exchangeConn = None
            except:
                print ("Error on closing sftp connection")
                # print sys.exc_info()
                # traceback.print_exc(file=sys.stdout)
                sys.stdout.flush()
        else:
            return "It is not opened!"

    def create_remote_link(self, fromfile, remotefile, stayalive=True):
        self.lock.acquire()
        assert self.isExchangeFileOpen
        buff = ""
        try:
            nomeu = "./" + os.path.basename(remotefile)
            if self.isnfs:
                # print "======",self.exchangeConn.getcwd()
                # print fromfile,nomeu
                o = self.exchangeConn.symlink(fromfile, nomeu)
                """
                while 1:
                        try:
                                filestat=self.exchangeConn.stat(nomeu)
                                #print "from frontend:", filestat
                                break
                        except:
                                #print "From frontend: Cannot stat "+str(nomeu)+" Trying again..."
                                timep.sleep(3)
                esegui = True
                while esegui:
                        out = self.connection.send_command_to_channel(self.channel,'stat '+str(os.path.join(self.get_remote_pwd(),os.path.basename(remotefile))),self.promptB)
                        #print "From submitter\n",out
                        outlines = out.split()

                        for word in outlines:
                                if not word.strip().startswith("stat:"):
                                        esegui = False
                                        break
                """
            else:
                buff += self.connection.send_command_to_channel(self.sftpsubcha, "symlink " + fromfile + " " + nomeu,
                                                                self.promptSFTP)
                # buff += self.connection.send_command_to_channel(self.channel,"symlink "+fromfile+" "+nomeu,self.promptB)
                """
                esegui = True
                while esegui:
                        out = self.connection.send_command_to_channel(self.channel,'stat '+str(nomeu),self.promptB)
                        print out
                        outlines = out.splitlines()
                        for linea in outlines:
                                if not linea.strip().startswith("stat:"):
                                        esegui = False
                                        break
                """
                # print buff
        except:
            print("Error on create a remote link")
            # print sys.exc_info()
            # traceback.print_exc()
        if not stayalive:
            buff += self.close_exchange_file()

        self.lock.release()
        return buff

    def create_remote_dir(self, remotedir, stayalive=True):
        self.lock.acquire()
        assert self.isExchangeFileOpen
        buff = ""
        try:
            if len(remotedir) >= 1 and remotedir[0] != "/" and not remotedir.startswith("./"):
                remotedir = "./" + remotedir

            buff += self.remove_remote_dir(remotedir)
            if self.isnfs:
                # print "mkdir -p ",remotedir
                self.exchangeConn.mkdir(remotedir)
                # buff += self.connection.send_command_to_channel(self.channel,"mkdir -p "+remotedir,self.promptB)
            else:
                buff += self.connection.send_command_to_channel(self.channel, "mkdir -p " + remotedir, self.promptB)
        except:
            print("Error on create remote directory", remotedir)
            # print sys.exc_info()
            # traceback.print_exc()
            # print "********************"
            # print buff
            # print "********************"

        if not stayalive:
            buff += self.close_exchange_file()

        self.lock.release()
        return buff

    def change_remote_dir(self, remotedir, stayalive=True):
        self.lock.acquire()
        assert self.isExchangeFileOpen
        buff = ""
        try:
            if len(remotedir) >= 1 and remotedir[0] != "/" and not remotedir.startswith("./"):
                remotedir = "./" + remotedir

            if self.isnfs:
                # print "cd ",remotedir
                self.exchangeConn.chdir(remotedir)
                # timep.sleep(10)
                # NOTE: we do not use yy. We really need to take it?
                yy = self.get_remote_listdir()
                buff += self.connection.send_command_to_channel(self.channel,
                                                                "cd " + str(os.path.abspath(self.get_remote_pwd())),
                                                                self.promptB)
            else:
                buff += self.connection.send_command_to_channel(self.sftpsubcha, "cd " + remotedir, self.promptSFTP)
                buff += self.connection.send_command_to_channel(self.channel, "cd " + remotedir, self.promptB)

                # print buff
        except:
            print("Error on change working remote directory", remotedir)
            # print sys.exc_info()
            # traceback.print_exc()
        if not stayalive:
            buff += self.close_exchange_file()

        self.lock.release()
        return buff

    def remove_remote_dir(self, remotedir, stayalive=True):
        self.lock.acquire()
        assert self.isExchangeFileOpen
        buff = ""
        try:
            remotedir = "./" + os.path.basename(os.path.normpath(remotedir))

            if self.isnfs:
                # print "rm -rf ",remotedir
                # self.exchangeConn.rmdir(remotedir)
                buff += self.connection.send_command_to_channel(self.channel, "rm -rf " + remotedir, self.promptB)
            else:
                buff += self.connection.send_command_to_channel(self.channel, "rm -rf " + remotedir, self.promptB)
        except:
            print("Error on remove remote directory maybe it does not exist")
            # print sys.exc_info()
            # traceback.print_exc()
        if not stayalive:
            buff += self.close_exchange_file()

        self.lock.release()
        return buff

    def create_remote_file(self, remotefile, stayalive=True):
        self.lock.acquire()
        assert self.isExchangeFileOpen
        buff = ""
        try:
            if self.isnfs:
                a = self.exchangeConn.open(remotefile, mode='w')
                a.close()
            else:
                buff += self.connection.send_command_to_channel(self.sftpsubcha, "!touch " + remotefile,
                                                                self.promptSFTP)
                buff += self.connection.send_command_to_channel(self.sftpsubcha, "put " + remotefile, self.promptSFTP)
                buff += self.connection.send_command_to_channel(self.sftpsubcha, "!rm " + remotefile, self.promptSFTP)
                # print buff
        except:
            print("Error on create remote file")
            # print sys.exc_info()
            # traceback.print_exc()
        if not stayalive:
            buff += self.close_exchange_file()
        self.lock.release()
        return buff

    def remove_remote_file(self, remotefile, stayalive=True):
        self.lock.acquire()
        assert self.isExchangeFileOpen
        buff = ""
        try:
            if self.isnfs:
                self.exchangeConn.remove(remotefile)
            else:
                buff += self.connection.send_command_to_channel(self.sftpsubcha, "rm -f" + remotefile, self.promptSFTP)
        except:
            print("Error on delete remote file", remotefile)
            # print sys.exc_info()
            # traceback.print_exc(file=sys.stdout)
        if not stayalive:
            buff += self.close_exchange_file()
        self.lock.release()
        return buff

    def copy_directory(self, localdir, remotedir, stayalive=True):
        self.lock.acquire()
        nomeu = "./" + os.path.basename(os.path.normpath(remotedir))
        buff = self.create_remote_dir(nomeu)
        buff += self.change_remote_dir(nomeu)
        tarro = tarfile.open(os.path.join(localdir, "transfert.tar.gz"), "w:gz")
        for root, subFolders, files in os.walk(localdir):
            for fileu in files:
                tarro.add(os.path.join(root, fileu), arcname="./" + fileu)
        tarro.close()

        buff += self.copy_local_file(os.path.join(localdir, "transfert.tar.gz"), "transfert.tar.gz")
        buff += self.connection.send_command_to_channel(self.channel, 'tar -zxf ' + str(
            os.path.join(self.get_remote_pwd(), 'transfert.tar.gz')), self.promptB)
        buff += self.remove_remote_file("transfert.tar.gz")
        os.remove(os.path.join(localdir, "transfert.tar.gz"))
        buff += self.change_remote_dir("..")

        if not stayalive:
            buff += self.close_exchange_file()
        self.lock.release()
        return buff

    def copy_local_file(self, localfile, remotefile, remote_path_asitis=False, stayalive=True, force_cumulative=False,
                        send_now=False):
        self.lock.acquire()
        aa = ""
        if send_now:
            aa = self.__real_copy_local_file(localfile, remotefile, stayalive=stayalive,
                                             remote_path_asitis=remote_path_asitis)
        elif os.path.basename(localfile)[-4:] in self.EXCEPTION_TO_CUMULATIVE and not force_cumulative:
            aa = self.__real_copy_local_file(localfile, remotefile, stayalive=stayalive,
                                             remote_path_asitis=remote_path_asitis)
        else:
            self.FILE_TO_COPY.append((localfile, remotefile))
            aa = "Appending transfert file request for " + localfile + "\n"
        self.lock.release()
        return aa

    def __real_copy_local_file(self, localfile, remotefile, stayalive=True, remote_path_asitis=False):
        assert self.isExchangeFileOpen
        buff = ""
        nowDir = ""
        try:
            if remote_path_asitis:
                nowDir = self.get_remote_pwd()
                self.change_remote_dir(os.path.dirname(remotefile))
            nomeu = "./" + os.path.basename(remotefile)
            if self.isnfs:
                # print "...........",self.exchangeConn.getcwd()
                o = self.exchangeConn.put(localfile, nomeu)
                while 1:
                    try:
                        filestat = self.exchangeConn.stat(nomeu)
                        # print "from frontend:", filestat
                        break
                    except:
                        # print "From frontend: Cannot stat "+str(nomeu)+" Trying again..."
                        timep.sleep(3)
                esegui = True
                while esegui:
                    out = self.connection.send_command_to_channel(self.channel, 'stat ' + str(
                        os.path.join(self.get_remote_pwd(), os.path.basename(remotefile))), self.promptB)
                    # print "From submitter\n",out
                    outlines = out.split()

                    for word in outlines:
                        if not word.strip().startswith("stat:"):  # "File:" is not ok for non english SO
                            esegui = False
                            break
            else:
                o = self.exchangeConn.put(localfile, nomeu)
                while 1:
                    try:
                        filestat = self.exchangeConn.stat(nomeu)
                        # print filestat
                        break
                    except:
                        print("Cannot stat " + str(nomeu) + " Trying again...")
                        timep.sleep(3)

                buff += self.connection.send_command_to_channel(self.sftpsubcha, "put " + nomeu, self.promptSFTP)
                buff += self.connection.send_command_to_channel(self.sftpsubcha, "!rm " + nomeu, self.promptSFTP)

                esegui = True
                while esegui:
                    out = self.connection.send_command_to_channel(self.channel, 'stat ' + str(nomeu), self.promptB)
                    # print out
                    outlines = out.splitlines()
                    for linea in outlines:
                        if not linea.strip().startswith("stat:"):
                            esegui = False
                            break

                            # print buff
        except:
            print("Error on putting local file into remote file")
            # print sys.exc_info()
            # traceback.print_exc()
        finally:
            if not stayalive:
                buff += self.close_exchange_file()

            if remote_path_asitis:
                self.change_remote_dir(nowDir)
            return buff

    def chmod_remote_file(self, remotefile, mode, stayalive=True):
        self.lock.acquire()
        assert self.isExchangeFileOpen
        buff = ""
        try:
            nomeu = "./" + os.path.basename(remotefile)
            if self.isnfs:
                self.exchangeConn.chmod(nomeu, mode)
            else:
                buff += self.connection.send_command_to_channel(self.sftpsubcha, "!chmod " + str(mode) + " " + nomeu,
                                                                self.promptSFTP)
                # print buff
        except:
            print("Error on getting remote file to local file")
            # print sys.exc_info()
            # traceback.print_exc()
        if not stayalive:
            buff += self.close_exchange_file()
        self.lock.release()
        return buff

    def get_remote_listdir(self, typef="*"):
        self.lock.acquire()
        # out = self.connection.send_command_to_channel(self.channel,'ls -1 '+str(typef)+' | sort -k1 -n',self.promptB)
        out = self.connection.send_command_to_channel(self.channel, 'echo $SHELL', self.promptB)
        shelltype = os.path.basename((out.splitlines()[-2]))
        if shelltype.upper() in ["TCSH"]:
            out = self.connection.send_command_to_channel(self.channel, 'foreach file (.)', "foreach? ")
            out = self.connection.send_command_to_channel(self.channel, 'ls -1 $file', "foreach? ")
            out = self.connection.send_command_to_channel(self.channel, 'end', self.promptB)
        elif shelltype.upper() in ["BASH", "SH"]:
            out = self.connection.send_command_to_channel(self.channel,
                                                          'for each_file in .; do echo "`ls -1 $each_file`"; done;',
                                                          self.promptB)
        listus = out.splitlines()
        listus = listus[1:-1]
        ab = []
        for li in listus:
            miax = li.strip()
            if typef == "*" and miax != "" and len(miax.split()) == 1:
                ab.append(miax)
            elif miax != "" and len(miax.split()) == 1 and miax.endswith(typef[1:]):
                ab.append(miax)
        self.lock.release()
        # print sorted(ab,self.__cmp_pdb)
        return sorted(ab, key=Grid.__cmp_pdb)

    @classmethod
    def __cmp_pdb(a, b):
        try:
            ai = a.split(".")[0]
            bi = b.split(".")[0]
            if len(ai.split("_")) > 1:
                ai = int(ai.split("_")[0])
            else:
                ai = int(ai)

            if len(bi.split("_")) > 1:
                bi = int(bi.split("_")[0])
            else:
                bi = int(bi)
        except:
            ai = a
            bi = b

        return (ai > bi) - (ai < bi)

    def get_remote_file(self, remotefile, localfile, stayalive=True, tryonetime=False, relaunch=False, command=None,
                        conditioEND=None, testEND=None, only_get_this=False, lenght_ext=4):
        # print "requiring",remotefile,"into",localfile
        self.lock.acquire()
        buff = ""
        if remotefile in self.FILE_RECEIVED:
            del self.FILE_RECEIVED[remotefile]
            self.lock.release()
            return "File: " + remotefile + " received!"
        else:
            lista_files = []
            while 1:
                try:
                    self.actualDirectory = self.get_remote_pwd()
                    lista_files = self.get_remote_listdir(typef="*" + remotefile[-1 * (lenght_ext):])
                    if len(lista_files) > 0 and remotefile in lista_files:
                        # print "remotefile is",remotefile
                        break
                    else:
                        # print "Nothing is ready yet, waiting 3 seconds..."
                        self.lock.release()
                        return False
                except:
                    print(sys.exc_info())
                    traceback.print_exc()
                    # RECONNECT
                    print("----- TRYING TO RECONNECT ------")
                    SystemUtility.remote_reconnection(self.actualDirectory)

            atleastone = False

            for cus in range(len(lista_files)):
                if cus > 300:
                    break
                fileu = lista_files[cus]
                base_fileu = os.path.basename(fileu)

                justcheck = True
                if base_fileu == remotefile:
                    justcheck = False
                elif not only_get_this:
                    justcheck = True
                else:
                    continue

                self.lock.release()
                # print "JUSTCHECK",justcheck,base_fileu,remotefile
                response = self.__real_get_remote_file(base_fileu, localfile, stayalive=stayalive,
                                                       tryonetime=tryonetime, relaunch=relaunch, command=command,
                                                       conditioEND=conditioEND, testEND=testEND, onlycheck=justcheck)
                self.lock.acquire()
                if isinstance(response, bool) and response:
                    out = self.connection.send_command_to_channel(self.channel, 'tar -rf ' + str(
                        os.path.join(self.get_remote_pwd(), '../receveid.tar')) + ' .././' + str(
                        os.path.basename(os.path.dirname(localfile))) + "/" + str(base_fileu), self.promptB)
                    self.FILE_RECEIVED[base_fileu] = 0
                    buff += self.remove_remote_file(fileu)
                    atleastone = True
                elif isinstance(response, str):
                    # print "IT SHOULD BE A JUSTCHECK FALSE"
                    buff += self.remove_remote_file(fileu)
                    # print "response: ",response
                    # else:
                    #       print "FILE: ",fileu,"give ",response
            if atleastone:
                out = self.connection.send_command_to_channel(self.channel, 'gzip -9 ../receveid.tar', self.promptB)
                buff += self.__real_get_remote_file("../receveid.tar.gz",
                                                    os.path.join(os.path.dirname(localfile), "../receveid.tar.gz"),
                                                    stayalive=stayalive)
                # buff += self.__real_get_remote_file("../receveid.tar.gz",os.path.join(os.path.dirname(localfile),"receveid.tar.gz"),stayalive=stayalive)
                buff += self.remove_remote_file("../receveid.tar.gz")
                # print "================================================",buff
                tarro = tarfile.open(os.path.join(os.path.dirname(localfile), "../receveid.tar.gz"), "r:gz")
                tarro.extractall(path=os.path.join(os.path.dirname(localfile), "../"))
                tarro.close()
                os.remove(os.path.join(os.path.dirname(localfile), "../receveid.tar.gz"))
                """
                tarro = tarfile.open(os.path.join(os.path.dirname(localfile),"receveid.tar.gz"), "r:gz")
                tarro.extractall(path=os.path.dirname(localfile))
                tarro.close()
                os.remove(os.path.join(os.path.dirname(localfile),"receveid.tar.gz"))
                """
            self.lock.release()
            # print "AAAAAAAAA",buff
            return buff

    def __real_get_remote_file(self, remotefile, localfile, stayalive=True, tryonetime=False, relaunch=False,
                               command=None, conditioEND=None, testEND=None, onlycheck=False):
        assert self.isExchangeFileOpen
        buff = ""
        try:
            if not remotefile.startswith("../"):
                nomeu = "./" + os.path.basename(remotefile)
            else:
                nomeu = remotefile

            if self.isnfs:
                while 1:
                    try:
                        self.lock.acquire()
                        filestat = self.exchangeConn.stat(nomeu)
                        self.lock.release()
                        # print "from frontend:", filestat
                        break
                    except:
                        # print "Cannot stat "+str(nomeu)+" Trying again..."
                        self.lock.release()
                        if tryonetime or onlycheck:
                            return False
                        # print "Vado in sleep A"
                        timep.sleep(3)
                esegui = True
                times = 0
                correct = False
                while esegui:
                    self.lock.acquire()
                    out = self.connection.send_command_to_channel(self.channel, 'stat ' + str(
                        os.path.join(self.get_remote_pwd(), nomeu)), self.promptB)
                    outlines = out.split()
                    self.lock.release()
                    for word in outlines:
                        if not word.strip().startswith("stat:"):
                            esegui = False
                            correct = True
                            break
                        if tryonetime:
                            esegui = False

                        if times == 10 and relaunch and command != None and not onlycheck:
                            # rilancia
                            self.lock.acquire()
                            out = self.connection.send_command_to_channel(self.channel, command, self.promptB)
                            self.lock.release()
                            timep.sleep(3)
                            # print "Vado in sleep B"
                            break
                    times += 1

                if onlycheck and not correct:
                    # print "RETURNING FALSE AT STAGE A"
                    return False
                elif tryonetime and not correct:
                    # print "RETURNING FALSE AT STAGE Q"
                    return False

                esegui = True
                if conditioEND != None:
                    conditioEND = conditioEND.replace("*", os.path.join(self.get_remote_pwd(), nomeu[:-4]))
                correct = False
                while esegui:
                    if conditioEND != None and testEND != None:
                        self.lock.acquire()
                        out = self.connection.send_command_to_channel(self.channel, conditioEND, self.promptB)
                        edere = out.splitlines()
                        # self.lock.notify()
                        self.lock.release()
                        # TEMPORARY
                        # print out
                        # print edere
                        # print "------------",edere[-2].strip(),str(testEND),edere[-2].strip() == str(testEND)
                        # print self.connection.send_command_to_channel(self.channel,"pwd\n",self.promptB)
                        # print edere[-2].strip(),str(testEND)
                        if edere[-2].strip() == str(testEND):
                            esegui = False
                            correct = True
                        else:
                            # print "Test falso non sono uguali vado in sleep"
                            if onlycheck:
                                # print "RETURNING FALSE AT STAGE MMMMMM"
                                return False
                            # print "Vado in sleep C"
                            timep.sleep(3)
                    else:
                        esegui = False
                        correct = True

                if onlycheck and not correct:
                    # print "RETURNING FALSE AT STAGE B"
                    return False

                if not onlycheck:
                    self.lock.acquire()
                    self.exchangeConn.get(nomeu, localfile)
                    self.lock.release()
                else:
                    # print "RETURNING TRUE AT STAGE MIAO"
                    return True
            else:
                times = 0
                esegui = True
                correct = False
                while esegui:
                    self.lock.acquire()
                    out = self.connection.send_command_to_channel(self.channel, 'stat ' + str(
                        os.path.join(self.get_remote_pwd(), nomeu)), self.promptB)
                    outlines = out.split()
                    self.lock.release()
                    for word in outlines:
                        if not word.strip().startswith("stat:"):
                            esegui = False
                            correct = True
                            break
                        if tryonetime:
                            esegui = False

                        if times == 10 and relaunch and command != None and not onlycheck:
                            # rilancia
                            self.lock.acquire()
                            out = self.connection.send_command_to_channel(self.channel, command, self.promptB)
                            self.lock.release()
                            # print "Vado in sleep D"
                            timep.sleep(3)
                            break
                    times += 1

                if onlycheck and not correct:
                    # print "RETURNING FALSE AT STAGE C"
                    return False
                elif tryonetime and not correct:
                    return False
                # print "Superato primo test"
                esegui = True
                if conditioEND != None:
                    conditioEND = conditioEND.replace("*", os.path.join(self.get_remote_pwd(), nomeu[:-4]))

                correct = False
                while esegui:
                    if conditioEND != None and testEND != None:
                        self.lock.acquire()
                        out = self.connection.send_command_to_channel(self.channel, conditioEND, self.promptB)
                        edere = out.splitlines()
                        self.lock.release()
                        # TEMPORARY
                        # print "uscito"
                        # print out
                        # print edere
                        if edere[-2].strip() == str(testEND):
                            esegui = False
                            correct = True
                        else:
                            if onlycheck:
                                # print "RETURNING FALSE AT STAGE OOOOO"
                                # print out
                                return False
                            # print "Vado in sleep E"
                            timep.sleep(3)
                    else:
                        esegui = False
                        correct = True

                if onlycheck and not correct:
                    # print "RETURNING FALSE AT STAGE D"
                    return False
                # print "Superato secondo test"
                if not onlycheck:
                    self.lock.acquire()
                    buff += self.connection.send_command_to_channel(self.sftpsubcha, "get " + nomeu, self.promptSFTP)
                    self.exchangeConn.get(os.path.basename(nomeu), localfile)
                    buff += self.connection.send_command_to_channel(self.sftpsubcha, "!rm " + os.path.basename(nomeu),
                                                                    self.promptSFTP)
                    self.lock.release()
                else:
                    # print "Ritorno True"
                    return True
        except:
            print("Error on getting remote file to local file")
            # print sys.exc_info()
            # traceback.print_exc()
            traceback.print_exc(file=sys.stdout)
        if not stayalive:
            self.lock.acquire()
            buff += self.close_exchange_file()
            self.lock.release()
        return buff

    def get_remote_pwd(self, stayalive=True):
        self.lock.acquire()
        # print "Ho il lock"
        assert self.isExchangeFileOpen
        # print "Supero l'assert"
        buff = ""
        current = ""
        try:
            if self.isnfs:
                current = self.exchangeConn.getcwd()
            else:
                buff += self.connection.send_command_to_channel(self.sftpsubcha, "pwd", self.promptSFTP)
                # print buff
                lion = buff.split()
                for pat in lion:
                    if pat[0] == "/":
                        current = pat
                        break
        except:
            print("Error on getting remote working directory")
            # print sys.exc_info()
            # traceback.print_exc()

        if not stayalive:
            buff += self.close_exchange_file()
        self.lock.release()
        return current

    def submitJob(self, job):
        pass

    def getStatus(self, jobid, nqueue):
        pass

    def getGridQueue(self, cluster="", queue=""):
        pass

    def isGridAlive(self):
        pass

    def removeJob(self, jobid, nqueue):
        pass

    def removeCluster(self, jobid):
        pass

    def submitJobs(self, job, nqueue):
        pass

    def setRequirements(self, requiString):
        pass

    def setMemory(self, memoryString):
        pass

    def setRank(self, rankString):
        pass

    def getCMD(self, jobid):
        pass


class SLURMManager(Grid):
    """Manager for SLURM jobs """

    def __init__(self, partition=''):
        super(SLURMManager, self).__init__()
        self.type_grid = "slurm"
        self.partition = partition

    def submitJob(self, job, isthelast=False, forcesubmit=False):
        return self.submitJobs(job, 1, is_array_job=False)

    def getStatus(self, jobid, nqueue):
        if jobid not in self.JOBS.keys() and jobid not in self.JOBS_DONE.keys():
            # print "JNP A version"
            return "JNP"  # Job Never Performed
        elif jobid not in self.JOBS.keys() and jobid in self.JOBS_DONE.keys():
            return "NRA"  # Not Registered Anymore
        elif jobid in self.JOBS.keys() and jobid not in self.JOBS_DONE.keys() and len(self.JOBS[jobid]) == 0:
            # print "JNP B version"
            return "JNP"  # Job Never Performed

        while True:
            try:
                cluster = ((self.JOBS[jobid])[nqueue])["gridCluster"]
                nq = ((self.JOBS[jobid])[-1])["queue"]
                break
            except:
                pass

        if nqueue < 0:
            return "NCV"  # Non Corresponding Values
        dizio = self.getGridQueue()
        if cluster not in dizio.keys():
            ntot = ((self.JOBS[jobid])[nqueue])["qtotal"]
            ((self.JOBS[jobid])[nqueue]) = {}
            completed = True
            for tr in self.JOBS[jobid]:
                if len(tr.keys()) > 0:
                    completed = False
                    break
            if completed and len(self.JOBS[jobid]) == ntot:
                self.JOBS_DONE[jobid] = copy.deepcopy(self.JOBS[jobid])
                del self.JOBS[jobid]
            return "NRA"  # Not Registered Anymore
        return ((self.JOBS[jobid])[nqueue], dizio[cluster])

    def getGridQueue(self, cluster="", queue=""):
        pass


    def isGridAlive(self):
        pass

    def removeJob(self, jobid, nqueue):
        pass

    def removeCluster(self, jobid):
        pass

    def submitJobs(self, job, nqueue, is_array_job=True):
        if not os.path.exists("./temp/"):
            os.makedirs("./temp/")

        if hasattr(self, "channel"):
            tarro = tarfile.open(os.path.join("./temp/", "transfert.tar.gz"), "w:gz")

            for from_fi, to_fi in self.FILE_TO_COPY:
                tarro.add(from_fi, arcname=to_fi)
            tarro.close()

            buff = self.copy_local_file(os.path.join("./temp/", "transfert.tar.gz"), "transfert.tar.gz")
            buff = self.connection.send_command_to_channel(self.channel, 'tar -zxf ' + str(
                os.path.join(self.get_remote_pwd(), 'transfert.tar.gz')), self.promptB)
            buff = self.remove_remote_file("transfert.tar.gz")
            os.remove(os.path.join("./temp/", "transfert.tar.gz"))
            self.FILE_TO_COPY = []

        while True:
            if len(threading.enumerate()) < 10:
                # print "Starting Thread Job..."
                new_t = SystemUtility.OutputThreading(self.__submitJobs, job, nqueue, is_array_job)
                new_t.start()
                # print "...started!"
                break
            else:
                print("Too many thread", len(threading.enumerate()))
                timep.sleep(3)
        return 0, nqueue

    def __submitJobs(self, job, nqueue, is_array_job):
        workDir = ""
        if hasattr(self, "channel"):
            workDir = self.get_remote_pwd()

        print("WORKING DIRECTORY:", workDir)
        if job.getName() in self.JOBS.keys():
            raise TypeError(
                'Another job is already registered with ' + job.getName() + " and it is impossible to register two jobs with the same name.")

        if not os.path.exists("./grid_jobs/"):
            os.makedirs("./grid_jobs/")

        all_cmds = []
        basecmd = """#! /bin/bash
# Template for SLURM array job for ARCIMBOLDO_BORGES and ARCIMBOLDO_LITE
#  Cambiar el valor de wd, el rango es -t 1-n, empieza en 1
# Una vez empezados los parametros de SLURM, no poner lineas en blanco o deja de leerlos.
"""

        text_partition = ""
        if self.partition != '':
            text_partition = "#SBATCH --partition=" + self.partition + "\n"
        if isinstance(job.initdir, str):
            basecmd += "#SBATCH --job-name=p" + str(job.getName()) + "\n" + text_partition + "#SBATCH --time=05:00:00\n"
            if is_array_job:
                basecmd += "#SBATCH --array=0-" + str(nqueue - 1) + "\n" + "#SBATCH --ntasks=1\n\n"
        elif isinstance(job.initdir, list):
            basecmd += "#SBATCH --job-name=p" + str(
                job.getName()) + "\n" + text_partition + "#SBATCH --time=05:00:00\n#SBATCH --ntasks=1\n\n"

        last_dire = ""
        glcount = 1

        for i in range(nqueue):
            cmd = basecmd
            commandline = ""
            if isinstance(job.initdir, str):
                if job.initdir != "":
                    commandline = "-D " + job.initdir + " -o /dev/null -e /dev/null"
                else:
                    commandline = "-o /dev/null -e /dev/null --requeue"
                if i > 0:
                    break
            elif isinstance(job.initdir, list):
                dire = ""
                nume = 0
                summa = 0
                for el in job.initdir:
                    dire, nume = el
                    summa += nume
                    if i < summa:
                        break
                # print "last_dire",last_dire,dire
                if last_dire != dire:
                    cmd += "#SBATCH --array=" + str(i) + "-" + str(summa - 1) + "\n\n"
                    commandline = "-D " + str(dire) + " -o /dev/null -e /dev/null --requeue"
                    last_dire = dire
                else:
                    continue

            # print "actual cmd"
            # print cmd
            cmd += "\nsrun "
            if not hasattr(self, "channel") and job.executable.endswith(".py"):
                cmd += '' + PATH_LOCAL_PYTHON_INTERPRETER + ' ' + job.executable
            elif job.executable.endswith(".py"):
                cmd += '' + PATH_REMOTE_PYTHON_INTERPRETER + ' ' + job.executable
            else:
                cmd += '' + job.executable

            if len(job.args) > 0:
                data = job.getArgs(True)
                if is_array_job:
                    # NOTE CM: changing the replacement environment variable to avoid errors
                    #data = data.replace("$(Process)", "$SLURM_ARRAY_TASK_ID")
                    data = data.replace("$(Process)", "${SLURM_ARRAY_TASK_ID}")
                else:
                    data = data.replace("$(Process)", "")

                cmd += " " + data

            if len(job.stdIn) > 0:
                data = job.getStdIn(True)
                if is_array_job:
                    #data = data.replace("$(Process)", "$SLURM_ARRAY_TASK_ID")
                    data = data.replace("$(Process)", "${SLURM_ARRAY_TASK_ID}")
                else:
                    data = data.replace("$(Process)", "")

                cmd += "<" + data

            if len(job.listFilesOut) > 1 and job.getOutput(True, "error") == job.getOutput(True, "output"):
                data = job.getOutput(True, "output")
                if is_array_job:
                    #data = data.replace("$(Process)", "$SLURM_ARRAY_TASK_ID")
                    data = data.replace("$(Process)", "${SLURM_ARRAY_TASK_ID}")
                else:
                    data = data.replace("$(Process)", "")

                cmd += ">&" + data
            elif len(job.listFilesOut) > 0:
                data = job.getOutput(True, "output")
                if is_array_job:
                    #data = data.replace("$(Process)", "$SLURM_ARRAY_TASK_ID")
                    data = data.replace("$(Process)", "${SLURM_ARRAY_TASK_ID}")
                else:
                    data = data.replace("$(Process)", "")

                cmd += ">" + data

            if len(job.listFilesOut) > 1 and job.getOutput(True, "error") != job.getOutput(True, "output"):
                data = job.getOutput(True, "error")
                if is_array_job:
                    #data = data.replace("$(Process)", "$SLURM_ARRAY_TASK_ID")
                    data = data.replace("$(Process)", "${SLURM_ARRAY_TASK_ID}")
                else:
                    data = data.replace("$(Process)", "")

                cmd += "2>" + data

            cmd += '\n'

            f = open("./grid_jobs/" + job.getName() + "_" + str(glcount) + ".sh", "w")
            f.write(cmd)
            f.close()
            all_cmds.append(("./grid_jobs/" + job.getName() + "_" + str(glcount) + ".sh", commandline))
            glcount += 1

        sent = 0
        tosend = nqueue
        self.JOBS[job.getName()] = []
        nproc = 0

        timep.sleep(60)  # This should avoid NFS delay errors in telemachus
        for cmde in all_cmds:
            cmd, commandline = cmde
            try:
                if not hasattr(self, "channel"):
                    nproc = tosend
                    while True:
                        p = subprocess.Popen(['sbatch'] + commandline.split() + [cmd], stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE)
                        print("sbatch " + commandline + " " + cmd.replace('"', '').strip())
                        out, err = p.communicate()
                        out = out.decode("ascii")
                        err = err.decode("ascii")
                        # scrivere nel Log l'eventuale err.
                        # print "======================OUTPUT===================="
                        print(out)
                        # print "======================OUTPUT===================="
                        # print "======================ERROR===================="
                        print(err)
                        # print "======================ERROR===================="
                        entered = False
                        liout = out.splitlines()
                        for lineacorr in liout:
                            try:
                                lir = int(lineacorr.strip().split()[3])
                                entered = True
                                break
                            except:
                                continue

                        if not entered:
                            print("SLURM not available sleeping 10 seconds...")
                            timep.sleep(10)
                            continue
                        nQ = sent
                        nCluster = int(lineacorr.strip().split()[3])
                        now = datetime.datetime.now()
                        time = now.strftime("%Y-%m-%d %H:%M")
                        # self.JOBS[job.getName()].append({"timeStart":time,"cmd":cmd,"job":job,"gridCluster":nCluster,"qtotal":nqueue,"queue":nQ})
                        tosend -= 1
                        sent += 1
                        nproc -= 1
                        print("job is: ", nCluster, nQ)
                        break
                else:
                    nproc = tosend
                    # print "Request for lock in submit",cmd
                    print("COPYING", cmd, "INTO", os.path.join(workDir, os.path.basename(cmd)))
                    self.copy_local_file(cmd, os.path.join(workDir, os.path.basename(cmd)), send_now=True,
                                         remote_path_asitis=True)
                    while True:
                        self.lock.acquire()
                        print("EXECUTING", 'sbatch ' + commandline + " " + os.path.join(workDir, os.path.basename(cmd)))
                        out = self.connection.send_command_to_channel(self.channel,
                                                                      'sbatch ' + commandline + " " + os.path.join(
                                                                          workDir, os.path.basename(cmd)), self.promptB)
                        self.lock.release()
                        # print "Lock free"
                        # scrivere nel Log l'eventuale err.
                        # print cmd
                        # print "======================OUTPUT===================="
                        # print out
                        # print "======================OUTPUT===================="
                        entered = False
                        liout = out.splitlines()
                        for lineacorr in liout:
                            try:
                                lir = int(lineacorr.strip().split()[3])
                                entered = True
                                break
                            except:
                                continue

                        if not entered:
                            print("SLURM not available sleeping 10 seconds...")
                            timep.sleep(10)
                            continue
                        nQ = sent
                        nCluster = int(lineacorr.strip().split()[3])
                        now = datetime.datetime.now()
                        time = now.strftime("%Y-%m-%d %H:%M")
                        # self.JOBS[job.getName()].append({"timeStart":time,"cmd":cmd,"job":job,"gridCluster":nCluster,"qtotal":nqueue, "queue":nQ})
                        tosend -= 1
                        sent += 1
                        nproc -= 1
                        print("job is: ", nCluster, nQ)
                        break
            except:
                print("Error on submitting a job...")
                print(sys.exc_info())
                traceback.print_exc(file=sys.stdout)
                timep.sleep(3)

    def setRequirements(self, requiString):
        pass

    def setMemory(self, memoryString):
        pass

    def setRank(self, rankString):
        pass

    def getCMD(self, jobid):
        try:
            return (self.JOBS[jobid])["cmd"]
        except:
            raise TypeError('The Job is not registered to the GridManager')


class MOABManager(Grid):
    """Manager for MOAB jobs """

    def __init__(self, partition='', fraction=1.0):
        super(MOABManager, self).__init__()
        self.type_grid = "moab"
        self.partition = partition
        self.fraction = fraction

    def submitJob(self, job, isthelast=False, forcesubmit=False):
        return self.submitJobs(job, 1, is_array_job=False)

    def getStatus(self, jobid, nqueue):
        if jobid not in self.JOBS.keys() and jobid not in self.JOBS_DONE.keys():
            # print "JNP A version"
            return "JNP"  # Job Never Performed
        elif jobid not in self.JOBS.keys() and jobid in self.JOBS_DONE.keys():
            return "NRA"  # Not Registered Anymore
        elif jobid in self.JOBS.keys() and jobid not in self.JOBS_DONE.keys() and len(self.JOBS[jobid]) == 0:
            # print "JNP B version"
            return "JNP"  # Job Never Performed

        while True:
            try:
                cluster = ((self.JOBS[jobid])[nqueue])["gridCluster"]
                nq = ((self.JOBS[jobid])[-1])["queue"]
                break
            except:
                pass

        if nqueue < 0:
            return "NCV"  # Non Corresponding Values
        dizio = self.getGridQueue()
        if cluster not in dizio.keys():
            ntot = ((self.JOBS[jobid])[nqueue])["qtotal"]
            ((self.JOBS[jobid])[nqueue]) = {}
            completed = True
            for tr in self.JOBS[jobid]:
                if len(tr.keys()) > 0:
                    completed = False
                    break
            if completed and len(self.JOBS[jobid]) == ntot:
                self.JOBS_DONE[jobid] = copy.deepcopy(self.JOBS[jobid])
                del self.JOBS[jobid]
            return "NRA"  # Not Registered Anymore
        return ((self.JOBS[jobid])[nqueue], dizio[cluster])

    def getGridQueue(self, cluster="", queue=""):
        pass

    def isGridAlive(self):
        pass

    def removeJob(self, jobid, nqueue):
        pass

    def removeCluster(self, jobid):
        pass

    def submitJobs(self, job, nqueue, is_array_job=True):
        if not os.path.exists("./temp/"):
            os.makedirs("./temp/")

        if hasattr(self, "channel"):
            tarro = tarfile.open(os.path.join("./temp/", "transfert.tar.gz"), "w:gz")

            for from_fi, to_fi in self.FILE_TO_COPY:
                tarro.add(from_fi, arcname=to_fi)
            tarro.close()

            buff = self.copy_local_file(os.path.join("./temp/", "transfert.tar.gz"), "transfert.tar.gz")
            buff = self.connection.send_command_to_channel(self.channel, 'tar -zxf ' + str(
                os.path.join(self.get_remote_pwd(), 'transfert.tar.gz')), self.promptB)
            buff = self.remove_remote_file("transfert.tar.gz")
            os.remove(os.path.join("./temp/", "transfert.tar.gz"))
            self.FILE_TO_COPY = []

        while True:
            if len(threading.enumerate()) < 10:
                # print "Starting Thread Job..."
                new_t = SystemUtility.OutputThreading(self.__submitJobs, job, nqueue, is_array_job)
                new_t.start()
                # print "...started!"
                break
            else:
                print("Too many thread", len(threading.enumerate()))
                timep.sleep(3)
        return 0, nqueue

    def __submitJobs(self, job, nqueue, is_array_job):
        workDir = ""
        if hasattr(self, "channel"):
            workDir = self.get_remote_pwd()

        print("WORKING DIRECTORY:", workDir)
        if job.getName() in self.JOBS.keys():
            raise TypeError(
                'Another job is already registered with ' + job.getName() + " and it is impossible to register two jobs with the same name.")

        if not os.path.exists("./grid_jobs/"):
            os.makedirs("./grid_jobs/")

        all_cmds = []
        basecmd = """
#! /bin/bash
# Template for SGE array job for ARCIMBOLDO-BORGES and ARCIMBOLDO
#  Cambiar el valor de wd, el rango es -t 1-n, empieza en 1
# Una vez empezados los parametros de SGE, no poner lineas en blanco o deja de leerlos.
"""

        if isinstance(job.initdir, str):
            basecmd += "#SBATCH --job-name=p" + str(
                job.getName()) + "\n" + "#SBATCH --partition=" + self.partition + "\n#SBATCH --time=05:00:00\n"
            if is_array_job:
                basecmd += "#SBATCH --array=0-" + str(nqueue - 1) + "\n" + "#SBATCH --ntasks=1\n\n"
        elif isinstance(job.initdir, list):
            basecmd += "#SBATCH --job-name=p" + str(
                job.getName()) + "\n" + "#SBATCH --partition=" + self.partition + "\n" + "#SBATCH --time=05:00:00\n#SBATCH --ntasks=1\n\n"

        last_dire = ""
        glcount = 1

        for i in range(nqueue):
            cmd = basecmd
            commandline = ""
            if isinstance(job.initdir, str):
                if job.initdir != "":
                    commandline = "-d " + job.initdir + " -o /dev/null -e /dev/null -l walltime=1:00:00:00"
                else:
                    commandline = "-o /dev/null -e /dev/null -r y -l walltime=1:00:00:00"
                if i > 0:
                    break
            elif isinstance(job.initdir, list):
                dire = ""
                nume = 0
                summa = 0
                for el in job.initdir:
                    dire, nume = el
                    summa += nume
                    if i < summa:
                        break
                # print "last_dire",last_dire,dire
                if last_dire != dire:
                    cmd += "#SBATCH --array=" + str(i) + "-" + str(summa - 1) + "\n\n"
                    commandline = "-d " + str(dire) + " -o /dev/null -e /dev/null -r y -l walltime=1:00:00:00"
                    last_dire = dire
                else:
                    continue

            # print "actual cmd"
            # print cmd
            if not hasattr(self, "channel") and job.executable.endswith(".py"):
                cmd += '' + PATH_LOCAL_PYTHON_INTERPRETER + ' ' + job.executable
            elif job.executable.endswith(".py"):
                cmd += '' + PATH_REMOTE_PYTHON_INTERPRETER + ' ' + job.executable
            else:
                cmd += '' + job.executable

            if len(job.args) > 0:
                data = job.getArgs(True)
                if is_array_job:
                    data = data.replace("$(Process)", "$SLURM_ARRAY_TASK_ID")
                else:
                    data = data.replace("$(Process)", "")

                cmd += " " + data

            if len(job.stdIn) > 0:
                data = job.getStdIn(True)
                if is_array_job:
                    data = data.replace("$(Process)", "$SLURM_ARRAY_TASK_ID")
                else:
                    data = data.replace("$(Process)", "")

                cmd += "<" + data

            if len(job.listFilesOut) > 1 and job.getOutput(True, "error") == job.getOutput(True, "output"):
                data = job.getOutput(True, "output")
                if is_array_job:
                    data = data.replace("$(Process)", "$SLURM_ARRAY_TASK_ID")
                else:
                    data = data.replace("$(Process)", "")

                cmd += ">&" + data
            elif len(job.listFilesOut) > 0:
                data = job.getOutput(True, "output")
                if is_array_job:
                    data = data.replace("$(Process)", "$SLURM_ARRAY_TASK_ID")
                else:
                    data = data.replace("$(Process)", "")

                cmd += ">" + data

            if len(job.listFilesOut) > 1 and job.getOutput(True, "error") != job.getOutput(True, "output"):
                data = job.getOutput(True, "error")
                if is_array_job:
                    data = data.replace("$(Process)", "$SLURM_ARRAY_TASK_ID")
                else:
                    data = data.replace("$(Process)", "")

                cmd += "2>" + data

            cmd += '\n'

            f = open("./grid_jobs/" + job.getName() + "_" + str(glcount) + ".sh", "w")
            f.write(cmd)
            f.close()
            all_cmds.append(("./grid_jobs/" + job.getName() + "_" + str(glcount) + ".sh", commandline))
            glcount += 1

        sent = 0
        tosend = nqueue
        self.JOBS[job.getName()] = []
        nproc = 0
        timep.sleep(60)  # This should avoid NFS delay errors in telemachus
        for cmde in all_cmds:
            cmd, commandline = cmde
            try:
                if not hasattr(self, "channel"):
                    if nproc == 0:
                        while True:
                            if float(self.fraction) == 1.0:
                                nproc = tosend
                                break

                            #(canI, nproc) = canSubmitJobs(getpass.getuser(), qname=self.qname, fraction=self.fraction)

                            if nproc >= 1:
                                if tosend < nproc:
                                    nproc = tosend
                                break
                            else:
                                timep.sleep(3)
                    while True:
                        p = subprocess.Popen(['msub'] + commandline.split() + [cmd], stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE)
                        # print cmd.replace('"','').strip()
                        out, err = p.communicate()
                        out = out.decode("ascii")
                        err = err.decode("ascii")
                        # scrivere nel Log l'eventuale err.
                        # print "======================OUTPUT===================="
                        print(out)
                        # print "======================OUTPUT===================="
                        # print "======================ERROR===================="
                        print(err)
                        # print "======================ERROR===================="
                        entered = False
                        liout = out.splitlines()
                        for lineacorr in liout:
                            try:
                                lir = int(lineacorr.strip())
                                entered = True
                                break
                            except:
                                continue

                        if not entered:
                            print("MOAB not available sleeping 10 seconds...")
                            timep.sleep(10)
                            continue
                        nQ = sent
                        nCluster = int(lineacorr.split(".")[0])
                        now = datetime.datetime.now()
                        time = now.strftime("%Y-%m-%d %H:%M")
                        # self.JOBS[job.getName()].append({"timeStart":time,"cmd":cmd,"job":job,"gridCluster":nCluster,"qtotal":nqueue,"queue":nQ})
                        tosend -= 1
                        sent += 1
                        nproc -= 1
                        print("job is: ", nCluster, nQ)
                        break
                else:
                    if nproc == 0:
                        while True:
                            if float(self.fraction) == 1.0:
                                nproc = tosend
                                break

                            self.lock.acquire()
                            out = self.connection.send_command_to_channel(self.channel,
                                                                          PATH_REMOTE_SGEPY + ' -qsub ' + str(
                                                                              self.qname) + " " + str(
                                                                              self.fraction) + "\n", self.promptB)
                            self.lock.release()

                            lis = out.splitlines()
                            lineacorr = ""
                            for lineacorr in lis:
                                if lineacorr.strip().startswith("Can"):
                                    break

                            canI = bool(lineacorr.split()[1])
                            nproc = int(lineacorr.split()[3])

                            # print "I can process:",canI,nproc

                            # nproc = 4 #TEMPORARY FOR TESTING

                            if nproc >= 1:
                                if tosend < nproc:
                                    nproc = tosend
                                break
                            else:
                                timep.sleep(3)
                    # print "Request for lock in submit",cmd
                    print("COPYING", cmd, "INTO", os.path.join(workDir, os.path.basename(cmd)))
                    self.copy_local_file(cmd, os.path.join(workDir, os.path.basename(cmd)), send_now=True,
                                         remote_path_asitis=True)
                    while True:
                        self.lock.acquire()
                        print("EXECUTING", 'msub ' + commandline + " " + os.path.join(workDir, os.path.basename(cmd)))
                        out = self.connection.send_command_to_channel(self.channel,
                                                                      'msub ' + commandline + " " + os.path.join(
                                                                          workDir, os.path.basename(cmd)), self.promptB)
                        self.lock.release()
                        # print "Lock free"
                        # scrivere nel Log l'eventuale err.
                        # print cmd
                        # print "======================OUTPUT===================="
                        # print out
                        # print "======================OUTPUT===================="
                        entered = False
                        liout = out.splitlines()
                        for lineacorr in liout:
                            try:
                                lir = int(lineacorr.strip())
                                entered = True
                                break
                            except:
                                continue

                        if not entered:
                            print("MOAB not available sleeping 10 seconds...")
                            timep.sleep(10)
                            continue
                        nQ = sent
                        nCluster = int(lineacorr.split(".")[0])
                        now = datetime.datetime.now()
                        time = now.strftime("%Y-%m-%d %H:%M")
                        # self.JOBS[job.getName()].append({"timeStart":time,"cmd":cmd,"job":job,"gridCluster":nCluster,"qtotal":nqueue, "queue":nQ})
                        tosend -= 1
                        sent += 1
                        nproc -= 1
                        print("job is: ", nCluster, nQ)
                        break
            except:
                print("Error on submitting a job...")
                print(sys.exc_info())
                traceback.print_exc(file=sys.stdout)
                timep.sleep(3)

    def setRequirements(self, requiString):
        pass

    def setMemory(self, memoryString):
        pass

    def setRank(self, rankString):
        pass

    def getCMD(self, jobid):
        try:
            return (self.JOBS[jobid])["cmd"]
        except:
            raise TypeError('The Job is not registered to the GridManager')


class TORQUEManager(Grid):
    """Manager for TORQUE jobs """

    def __init__(self, qname='all.q', cores_per_node=4, parallel_jobs=100, maui=False):
        super(TORQUEManager, self).__init__()
        self.type_grid = "torque"
        self.qname = str(qname)
        self.maui = maui
        # while parallel_jobs%cores_per_node != 0 and parallel_jobs<1000:
        #    parallel_jobs += 1

        if parallel_jobs > 1000:
            print("ERROR: The minimum number of parallel jobs ", parallel_jobs, "that can be efficiently grouped in",
                  cores_per_node,
                  "cores per node is greater than 1000. for security reason is not allowed to send all those jobs in parallel! Contact the support of the ARCIMBOLDO Team, please.")
            sys.exit(1)

        self.cores_per_node = cores_per_node
        self.parallel_jobs = parallel_jobs
        self.node_maximum = int(self.parallel_jobs / self.cores_per_node)

    def submitJob(self, job, isthelast=False, forcesubmit=False):
        if not os.path.exists("./temp/"):
            os.makedirs("./temp/")

        if hasattr(self, "channel"):
            if not forcesubmit and len(self.FILE_TO_COPY) > 0 and (
                    len(self.FILE_TO_COPY) >= self.CUMULATIVE_TRANSFERING or isthelast):
                tarro = tarfile.open(os.path.join("./temp/", "transfert.tar.gz"), "w:gz")

                for from_fi, to_fi in self.FILE_TO_COPY:
                    tarro.add(from_fi, arcname=to_fi)
                tarro.close()

                buff = self.copy_local_file(os.path.join("./temp/", "transfert.tar.gz"), "transfert.tar.gz")
                buff = self.connection.send_command_to_channel(self.channel, 'tar -zxf ' + str(
                    os.path.join(self.get_remote_pwd(), 'transfert.tar.gz')), self.promptB)
                buff = self.remove_remote_file("transfert.tar.gz")
                os.remove(os.path.join("./temp/", "transfert.tar.gz"))
                self.FILE_TO_COPY = []
                self.JOB_TO_SUBMIT.append(job)
                for ju in self.JOB_TO_SUBMIT:
                    a = self.submitJob(ju, forcesubmit=True)
                self.JOB_TO_SUBMIT = []
                return 0, 0
            elif not forcesubmit and len(self.FILE_TO_COPY) > 0:
                self.JOB_TO_SUBMIT.append(job)
                return 0, 0

        while True:
            if len(threading.enumerate()) < 10:
                # print "Starting Thread Job..."
                new_t = SystemUtility.OutputThreading(self.__submitJob, job)
                new_t.start()
                # print "...started!"
                break
            else:
                print("Too many thread", len(threading.enumerate()))
                timep.sleep(3)
        return 0, 1

    def __submitJob(self, job):
        workDir = ""
        nqueue = 1
        if hasattr(self, "channel"):
            workDir = self.get_remote_pwd()

        print("WORKING DIRECTORY:", workDir)
        if job.getName() in self.JOBS.keys():
            raise TypeError(
                'Another job is already registered with ' + job.getName() + " and it is impossible to register two jobs with the same name.")

        if not os.path.exists("./grid_jobs/"):
            os.makedirs("./grid_jobs/")

        all_cmds = []
        basecmd = """
#! /bin/bash
# Template for TORQUE array job for ARCIMBOLDO-BORGES and ARCIMBOLDO
#  Cambiar el valor de wd, el rango es -t 1-n, empieza en 1
# Una vez empezados los parametros de SGE, no poner lineas en blanco o deja de leerlos.
"""
        if isinstance(job.initdir, str):
            if job.initdir != "":
                basecmd += "#PBS -d " + str(job.initdir) + "\n" + "#PBS -N p" + str(
                    job.getName()) + "\n" + "#PBS -q " + self.qname + "\n" + "#PBS -r y\n#PBS -o output.out\n#PBS -e error.err\n"
            else:
                basecmd += "#PBS -N p" + str(
                    job.getName()) + "\n" + "#PBS -q " + self.qname + "\n" + "#PBS -r y\n#PBS -o output.out\n#PBS -e error.err\n"

        elif isinstance(job.initdir, list):
            basecmd += "#PBS -N p" + str(
                job.getName()) + "\n" + "#PBS -q " + self.qname + "\n" + "#PBS -r y\n#PBS -o output.out\n#PBS -e error.err\n"

        basecmd += "#PBS -l nodes=1:ppn=1\n\n"
        last_dire = ""
        glcount = 1
        for i in range(nqueue):
            cmd = basecmd
            # print "actual cmd"
            # print cmd
            if not hasattr(self, "channel") and job.executable.endswith(".py"):
                cmd += ' ' + PATH_LOCAL_PYTHON_INTERPRETER + ' ' + job.executable
            elif job.executable.endswith(".py"):
                cmd += ' ' + PATH_REMOTE_PYTHON_INTERPRETER + ' ' + job.executable
            else:
                cmd += ' ' + job.executable

            if len(job.args) > 0:
                cmd += " " + job.getArgs(False)

            if len(job.stdIn) > 0:
                cmd += "<" + job.getStdIn(False)

            if len(job.listFilesOut) > 1 and job.getOutput(False, "error") == job.getOutput(False, "output"):
                data = job.getOutput(False, "output")
                cmd += ">&" + data
            elif len(job.listFilesOut) > 0:
                data = job.getOutput(False, "output")
                cmd += ">" + data

            if len(job.listFilesOut) > 1 and job.getOutput(False, "error") != job.getOutput(False, "output"):
                data = job.getOutput(False, "error")
                cmd += "2>" + data

            cmd += '\n'

            f = open("./grid_jobs/" + job.getName() + "_" + str(glcount) + ".sh", "w")
            f.write(cmd)
            f.close()
            all_cmds.append("./grid_jobs/" + job.getName() + "_" + str(glcount) + ".sh")
            glcount += 1

        sent = 0
        tosend = nqueue
        self.JOBS[job.getName()] = []
        nproc = 0

        timep.sleep(60)  # This should avoid NFS delay errors in telemachus
        for cmd in all_cmds:
            try:
                if not hasattr(self, "channel"):
                    while True:
                        p = subprocess.Popen(['qsub', cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        # print cmd.replace('"','').strip()
                        out, err = p.communicate()
                        out = out.decode("ascii")
                        err = err.decode("ascii")
                        # scrivere nel Log l'eventuale err.
                        # print "======================OUTPUT===================="
                        # print out
                        # print "======================OUTPUT===================="
                        # print "======================ERROR===================="
                        # print err
                        # print "======================ERROR===================="
                        entered = False
                        liout = out.splitlines()
                        for lineacorr in liout:
                            if lineacorr.strip().split(".")[0].isdigit():
                                entered = True
                                break
                        if not entered:
                            print("TORQUE not available sleeping 10 seconds...")
                            timep.sleep(10)
                            continue
                        nQ = sent
                        nCluster = int(lineacorr.strip().split(".")[0])
                        now = datetime.datetime.now()
                        time = now.strftime("%Y-%m-%d %H:%M")
                        # self.JOBS[job.getName()].append({"timeStart":time,"cmd":cmd,"job":job,"gridCluster":nCluster,"qtotal":nqueue,"queue":nQ})
                        tosend -= 1
                        sent += 1
                        nproc -= 1
                        print("job is: ", nCluster, nQ)
                        break
                else:
                    # print "Request for lock in submit",cmd
                    print("COPYING", cmd, "INTO", os.path.join(workDir, os.path.basename(cmd)))
                    self.copy_local_file(cmd, os.path.join(workDir, os.path.basename(cmd)), send_now=True,
                                         remote_path_asitis=True)
                    while True:
                        self.lock.acquire()
                        print("EXECUTING", 'nohup qsub ' + os.path.join(workDir,
                                                                        os.path.basename(cmd) + " | tee nohup.out"))
                        out = self.connection.send_command_to_channel(self.channel,
                                                                      'nohup qsub ' + os.path.join(workDir,
                                                                                                   os.path.basename(
                                                                                                       cmd)) + " | tee nohup.out",
                                                                      self.promptB)
                        self.lock.release()
                        # print "Lock free"
                        # scrivere nel Log l'eventuale err.
                        # print cmd
                        # print "======================OUTPUT===================="
                        # print out
                        # print "======================OUTPUT===================="
                        entered = False
                        liout = out.splitlines()
                        for lineacorr in liout:
                            if lineacorr.strip().split(".")[0].isdigit():
                                entered = True
                                break
                        if not entered:
                            print("TORQUE not available sleeping 10 seconds...")
                            timep.sleep(10)
                            continue
                        nQ = sent
                        nCluster = int(lineacorr.strip().split(".")[0])
                        now = datetime.datetime.now()
                        time = now.strftime("%Y-%m-%d %H:%M")
                        # self.JOBS[job.getName()].append({"timeStart":time,"cmd":cmd,"job":job,"gridCluster":nCluster,"qtotal":nqueue, "queue":nQ})
                        tosend -= 1
                        sent += 1
                        nproc -= 1
                        print("job is: ", nCluster, nQ)
                        break
            except:
                print("Error on submitting a job...")
                print(sys.exc_info())
                traceback.print_exc(file=sys.stdout)
                timep.sleep(3)

    def getStatus(self, jobid, nqueue):
        if jobid not in self.JOBS.keys() and jobid not in self.JOBS_DONE.keys():
            # print "JNP A version"
            return "JNP"  # Job Never Performed
        elif jobid not in self.JOBS.keys() and jobid in self.JOBS_DONE.keys():
            return "NRA"  # Not Registered Anymore
        elif jobid in self.JOBS.keys() and jobid not in self.JOBS_DONE.keys() and len(self.JOBS[jobid]) == 0:
            # print "JNP B version"
            return "JNP"  # Job Never Performed

        while True:
            try:
                cluster = ((self.JOBS[jobid])[nqueue])["gridCluster"]
                nq = ((self.JOBS[jobid])[-1])["queue"]
                break
            except:
                pass

        if nqueue < 0:
            return "NCV"  # Non Corresponding Values
        dizio = self.getGridQueue()
        if cluster not in dizio.keys():
            ntot = ((self.JOBS[jobid])[nqueue])["qtotal"]
            ((self.JOBS[jobid])[nqueue]) = {}
            completed = True
            for tr in self.JOBS[jobid]:
                if len(tr.keys()) > 0:
                    completed = False
                    break
            if completed and len(self.JOBS[jobid]) == ntot:
                self.JOBS_DONE[jobid] = copy.deepcopy(self.JOBS[jobid])
                del self.JOBS[jobid]
            return "NRA"  # Not Registered Anymore
        return ((self.JOBS[jobid])[nqueue], dizio[cluster])

    def getGridQueue(self, cluster="", queue=""):
        pass

    def isGridAlive(self):
        pass

    def removeJob(self, jobid, nqueue):
        pass

    def removeCluster(self, jobid):
        pass

    def submitJobs(self, job, nqueue):
        if not os.path.exists("./temp/"):
            os.makedirs("./temp/")

        if hasattr(self, "channel"):
            tarro = tarfile.open(os.path.join("./temp/", "transfert.tar.gz"), "w:gz")

            for from_fi, to_fi in self.FILE_TO_COPY:
                tarro.add(from_fi, arcname=to_fi)
            tarro.close()

            buff = self.copy_local_file(os.path.join("./temp/", "transfert.tar.gz"), "transfert.tar.gz")
            buff = self.connection.send_command_to_channel(self.channel, 'tar -zxf ' + str(
                os.path.join(self.get_remote_pwd(), 'transfert.tar.gz')), self.promptB)
            buff = self.remove_remote_file("transfert.tar.gz")
            os.remove(os.path.join("./temp/", "transfert.tar.gz"))
            self.FILE_TO_COPY = []

        while True:
            if len(threading.enumerate()) < 10:
                # print "Starting Thread Job..."
                new_t = SystemUtility.OutputThreading(self.__submitJobs, job, nqueue)
                new_t.start()
                # print "...started!"
                break
            else:
                print("Too many thread", len(threading.enumerate()))
                timep.sleep(3)
        return 0, nqueue

    def __submitJobs(self, job, nqueue):
        workDir = ""
        if hasattr(self, "channel"):
            workDir = self.get_remote_pwd()

        print("WORKING DIRECTORY:", workDir)
        if job.getName() in self.JOBS.keys():
            raise TypeError(
                'Another job is already registered with ' + job.getName() + " and it is impossible to register two jobs with the same name.")

        if not os.path.exists("./grid_jobs/"):
            os.makedirs("./grid_jobs/")

        all_cmds = []
        basecmd = """
#! /bin/bash
# Template for TORQUE array job for ARCIMBOLDO-BORGES and ARCIMBOLDO
#  Cambiar el valor de wd, el rango es -t 1-n, empieza en 1
# Una vez empezados los parametros de SGE, no poner lineas en blanco o deja de leerlos.
"""
        if isinstance(job.initdir, str):
            if job.initdir != "":
                basecmd += "#PBS -d " + str(job.initdir) + "\n" + "#PBS -N p" + str(
                    job.getName()) + "\n" + "#PBS -q " + self.qname + "\n" + "#PBS -r y\n#PBS -o output.out\n#PBS -e error.err\n"
            else:
                basecmd += "#PBS -N p" + str(
                    job.getName()) + "\n" + "#PBS -q " + self.qname + "\n" + "#PBS -r y\n#PBS -o output.out\n#PBS -e error.err\n"

        elif isinstance(job.initdir, list):
            basecmd += "#PBS -N p" + str(
                job.getName()) + "\n" + "#PBS -q " + self.qname + "\n" + "#PBS -r y\n#PBS -o output.out\n#PBS -e error.err\n"

        last_dire = ""
        cmd_last = ""
        glcount = 1
        for i in range(nqueue):
            cmd = basecmd
            if isinstance(job.initdir, str):
                if i > 0:
                    break
            elif isinstance(job.initdir, list):
                dire = ""
                nume = 0
                summa = 0
                for el in job.initdir:
                    dire, nume = el
                    summa += nume
                    if i < summa:
                        break
                # print "====================nqueue:",nqueue,"summa",summa,"dire",dire
                if last_dire != dire:
                    if (summa - i) == self.parallel_jobs:
                        cmd += "#PBS -d " + str(dire) + "\n" + "#PBS -t " + str(i + 1) + "-" + str(summa) + "\n"
                        if not self.maui:
                            cmd += "#PBS -l nodes=" + str(int(self.node_maximum)) + ":ppn=" + str(
                                int(self.cores_per_node)) + "\n\n"
                    else:
                        todo = int((summa - i) / self.cores_per_node)
                        todo_rest = (summa - i) % self.cores_per_node
                        if todo == 0:
                            cmd += "#PBS -d " + str(dire) + "\n" + "#PBS -t " + str(i + 1) + "-" + str(summa) + "\n"
                            if not self.maui:
                                cmd += "#PBS -l nodes=" + str(1) + ":ppn=" + str(summa) + "\n\n"
                        else:
                            if todo_rest == 0:
                                cmd += "#PBS -d " + str(dire) + "\n" + "#PBS -t " + str(i + 1) + "-" + str(summa) + "\n"
                                if not self.maui:
                                    cmd += "#PBS -l nodes=" + str(todo) + ":ppn=" + str(
                                        int(self.cores_per_node)) + "\n\n"
                            else:
                                cmd_last = cmd
                                cmd += "#PBS -d " + str(dire) + "\n" + "#PBS -t " + str(i + 1) + "-" + str(
                                    i + (todo * self.cores_per_node) + todo_rest) + "\n"
                                if not self.maui:
                                    cmd += "#PBS -l nodes=" + str(todo) + ":ppn=" + str(
                                        int(self.cores_per_node)) + "\n\n"
                                    # cmd += "#PBS -d "+str(dire)+"\n"+"#PBS -t "+str((todo*self.cores_per_node)+1)+"-"+str((todo*self.cores_per_node)+1+todo_rest)+"\n"
                                    # if not self.maui:
                                    #   cmd_last += "#PBS -l nodes="+str(1)+":ppn="+str(todo_rest)+"\n\n"

                    last_dire = dire
                else:
                    continue

            cmd_sec = "i=$(($PBS_ARRAYID - 1))\n"
            # print "actual cmd"
            # print cmd
            if not hasattr(self, "channel") and job.executable.endswith(".py"):
                cmd_sec += '' + PATH_LOCAL_PYTHON_INTERPRETER + ' ' + job.executable
            elif job.executable.endswith(".py"):
                cmd_sec += '' + PATH_REMOTE_PYTHON_INTERPRETER + ' ' + job.executable
            else:
                cmd_sec += '' + job.executable

            if len(job.args) > 0:
                data = job.getArgs(True)
                data = data.replace("$(Process)", "$i")

                cmd_sec += " " + data

            if len(job.stdIn) > 0:
                data = job.getStdIn(True)
                data = data.replace("$(Process)", "$i")

                cmd_sec += "<" + data

            if len(job.listFilesOut) > 1 and job.getOutput(True, "error") == job.getOutput(True, "output"):
                data = job.getOutput(True, "output")
                data = data.replace("$(Process)", "$i")

                cmd_sec += ">&" + data
            elif len(job.listFilesOut) > 0:
                data = job.getOutput(True, "output")
                data = data.replace("$(Process)", "$i")

                cmd_sec += ">" + data

            if len(job.listFilesOut) > 1 and job.getOutput(True, "error") != job.getOutput(True, "output"):
                data = job.getOutput(True, "error")
                data = data.replace("$(Process)", "$i")

                cmd_sec += "2>" + data

            cmd_sec += '\n'

            cmd += cmd_sec
            f = open("./grid_jobs/" + job.getName() + "_" + str(glcount) + ".sh", "w")
            f.write(cmd)
            f.close()
            all_cmds.append("./grid_jobs/" + job.getName() + "_" + str(glcount) + ".sh")
            glcount += 1

            if len(cmd_last) > 0:
                cmd_last += cmd_sec
                f = open("./grid_jobs/" + job.getName() + "_" + str(glcount) + ".sh", "w")
                f.write(cmd_last)
                f.close()
                all_cmds.append("./grid_jobs/" + job.getName() + "_" + str(glcount) + ".sh")
                glcount += 1

        sent = 0
        tosend = nqueue
        self.JOBS[job.getName()] = []
        nproc = 0

        timep.sleep(60)  # This should avoid NFS delay errors in telemachus
        for cmd in all_cmds:
            try:
                if not hasattr(self, "channel"):
                    while True:
                        p = subprocess.Popen(['qsub', cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        # print cmd.replace('"','').strip()
                        out, err = p.communicate()
                        out = out.decode("ascii")
                        err = err.decode("ascii")
                        # scrivere nel Log l'eventuale err.
                        # print "======================OUTPUT===================="
                        # print out
                        # print "======================OUTPUT===================="
                        # print "======================ERROR===================="
                        # print err
                        # print "======================ERROR===================="
                        entered = False
                        liout = out.splitlines()
                        for lineacorr in liout:
                            if lineacorr.strip().split(".")[0].split("[]")[0].isdigit():
                                entered = True
                                break
                        if not entered:
                            print("TORQUE not available sleeping 10 seconds...")
                            timep.sleep(10)
                            continue
                        nQ = sent
                        nCluster = int(lineacorr.strip().split(".")[0].split("[]")[0])
                        now = datetime.datetime.now()
                        time = now.strftime("%Y-%m-%d %H:%M")
                        # self.JOBS[job.getName()].append({"timeStart":time,"cmd":cmd,"job":job,"gridCluster":nCluster,"qtotal":nqueue,"queue":nQ})
                        tosend -= 1
                        sent += 1
                        nproc -= 1
                        print("job is: ", nCluster, nQ)
                        break
                else:
                    # print "Request for lock in submit",cmd
                    print("COPYING", cmd, "INTO", os.path.join(workDir, os.path.basename(cmd)))
                    self.copy_local_file(cmd, os.path.join(workDir, os.path.basename(cmd)), send_now=True,
                                         remote_path_asitis=True)
                    while True:
                        self.lock.acquire()
                        print("EXECUTING", 'nohup qsub ' + os.path.join(workDir,
                                                                        os.path.basename(cmd) + " | tee nohup.out"))
                        out = self.connection.send_command_to_channel(self.channel,
                                                                      'nohup qsub ' + os.path.join(workDir,
                                                                                                   os.path.basename(
                                                                                                       cmd)) + " | tee nohup.out",
                                                                      self.promptB)
                        self.lock.release()
                        # print "Lock free"
                        # scrivere nel Log l'eventuale err.
                        # print cmd
                        # print "======================OUTPUT===================="
                        # print out
                        # print "======================OUTPUT===================="
                        entered = False
                        liout = out.splitlines()
                        # print "=======",liout
                        for lineacorr in liout:
                            if lineacorr.strip().split(".")[0].split("[]")[0].isdigit():
                                entered = True
                                break
                        if not entered:
                            print("TORQUE not available sleeping 10 seconds...")
                            timep.sleep(10)
                            continue
                        nQ = sent
                        nCluster = int(lineacorr.strip().split(".")[0].split("[]")[0])
                        now = datetime.datetime.now()
                        time = now.strftime("%Y-%m-%d %H:%M")
                        # self.JOBS[job.getName()].append({"timeStart":time,"cmd":cmd,"job":job,"gridCluster":nCluster,"qtotal":nqueue, "queue":nQ})
                        tosend -= 1
                        sent += 1
                        nproc -= 1
                        print("job is: ", nCluster, nQ)
                        break
            except:
                print("Error on submitting a job...")
                print(sys.exc_info())
                traceback.print_exc(file=sys.stdout)
                timep.sleep(3)

    def setRequirements(self, requiString):
        pass

    def setMemory(self, memoryString):
        pass

    def setRank(self, rankString):
        pass

    def getCMD(self, jobid):
        try:
            return (self.JOBS[jobid])["cmd"]
        except:
            raise TypeError('The Job is not registered to the GridManager')


class SGEManager(Grid):
    """Manager for Sun Grid Engine jobs """

    def __init__(self, qname='all.q', fraction=1.0):
        super(SGEManager, self).__init__()
        self.type_grid = "sge"
        self.qname = qname
        self.fraction = fraction

    def submitJob(self, job, isthelast=False, forcesubmit=False):
        if not os.path.exists("./temp/"):
            os.makedirs("./temp/")

        if hasattr(self, "channel"):
            if not forcesubmit and len(self.FILE_TO_COPY) > 0 and (
                    len(self.FILE_TO_COPY) >= self.CUMULATIVE_TRANSFERING or isthelast):
                tarro = tarfile.open(os.path.join("./temp/", "transfert.tar.gz"), "w:gz")

                for from_fi, to_fi in self.FILE_TO_COPY:
                    tarro.add(from_fi, arcname=to_fi)
                tarro.close()

                buff = self.copy_local_file(os.path.join("./temp/", "transfert.tar.gz"), "transfert.tar.gz")
                buff = self.connection.send_command_to_channel(self.channel, 'tar -zxf ' + str(
                    os.path.join(self.get_remote_pwd(), 'transfert.tar.gz')), self.promptB)
                buff = self.remove_remote_file("transfert.tar.gz")
                os.remove(os.path.join("./temp/", "transfert.tar.gz"))
                self.FILE_TO_COPY = []
                self.JOB_TO_SUBMIT.append(job)
                for ju in self.JOB_TO_SUBMIT:
                    a = self.submitJob(ju, forcesubmit=True)
                self.JOB_TO_SUBMIT = []
                return 0, 0
            elif not forcesubmit and len(self.FILE_TO_COPY) > 0:
                self.JOB_TO_SUBMIT.append(job)
                return 0, 0

        while True:
            if len(threading.enumerate()) < 10:
                # print "Starting Thread Job..."
                new_t = SystemUtility.OutputThreading(self.__submitJob, job)
                new_t.start()
                # print "...started!"
                break
            else:
                print("Too many thread", len(threading.enumerate()))
                timep.sleep(3)
        return 0, 1

    def __submitJob(self, job):
        if job.getName() in self.JOBS.keys():
            raise TypeError(
                'Another job is already registered with ' + job.getName() + " and it is impossible to register two jobs with the same name.")

        if not os.path.exists("./grid_jobs/"):
            os.makedirs("./grid_jobs/")

        cmd = ""
        if job.initdir != "":
            cmd += "qsub -q " + self.qname + " -wd " + str(job.initdir) + " -r y -o /dev/null -e /dev/null -N p" + str(
                job.getName()) + " -b y"
        else:
            cmd += "qsub -q " + self.qname + " -cwd -o /dev/null -r y -e /dev/null -N p" + str(job.getName()) + " -b y"

        if not hasattr(self, "channel") and job.executable.endswith(".py"):
            cmd += ' "' + PATH_LOCAL_PYTHON_INTERPRETER + ' ' + job.executable
        elif job.executable.endswith(".py"):
            cmd += ' "' + PATH_REMOTE_PYTHON_INTERPRETER + ' ' + job.executable
        else:
            cmd += ' "' + job.executable

        if len(job.args) > 0:
            cmd += " " + job.getArgs(False)

        if len(job.stdIn) > 0:
            cmd += "<" + job.getStdIn(False)

        if len(job.listFilesOut) > 1 and job.getOutput(False, "error") == job.getOutput(False, "output"):
            data = job.getOutput(False, "output")
            cmd += ">&" + data
        elif len(job.listFilesOut) > 0:
            data = job.getOutput(False, "output")
            cmd += ">" + data

        if len(job.listFilesOut) > 1 and job.getOutput(False, "error") != job.getOutput(False, "output"):
            data = job.getOutput(False, "error")
            cmd += "2>" + data

        cmd += '"\n'

        self.JOBS[job.getName()] = []
        nproc = 0
        timep.sleep(60)  # This should avoid NFS delay errors in telemachus
        try:
            if not hasattr(self, "channel"):
                if nproc == 0:
                    while True:
                        if float(self.fraction) == 1.0:
                            nproc = 1
                            break

                        # TODO: Get the output and process it ti extract nproc
                        (canI, nproc) = self.canSubmitJobs(getpass.getuser(), qname=self.qname, fraction=self.fraction)
                        # nproc = 4 #TEMPORARY FOR TESTING

                        if nproc >= 1:
                            nproc = 1
                            break
                        else:
                            timep.sleep(3)
                while True:
                    p = subprocess.Popen(cmd.replace('"', '').strip().split(), stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)
                    out, err = p.communicate()
                    out = out.decode("ascii")
                    err = err.decode("ascii")
                    # scrivere nel Log l'eventuale err.
                    # print "======================OUTPUT===================="
                    print(out)
                    # print "======================OUTPUT===================="
                    # print "======================ERROR===================="
                    print(err)
                    # print "======================ERROR===================="

                    liout = out.splitlines()
                    entered = False
                    for lineacorr in liout:
                        if lineacorr.strip().startswith("Your"):
                            entered = True
                            break
                    if not entered:
                        print("SGE not available sleeping 10 seconds...")
                        timep.sleep(10)
                        continue

                    nQ = 0
                    #nCluster = int(lineacorr.split()[2])
                    nCluster = lineacorr
                    now = datetime.datetime.now()
                    time = now.strftime("%Y-%m-%d %H:%M")
                    # self.JOBS[job.getName()].append({"timeStart":time,"cmd":cmd,"job":job,"gridCluster":nCluster,"qtotal":1,"queue":nQ})
                    nproc -= 1
                    print("job is: ", nCluster, nQ)
                    break
            else:
                if nproc == 0:
                    while True:
                        if float(self.fraction) == 1.0:
                            nproc = 1
                            break
                        self.lock.acquire()
                        out = self.connection.send_command_to_channel(self.channel, PATH_REMOTE_SGEPY + ' -qsub ' + str(
                            self.qname) + " " + str(self.fraction) + "\n", self.promptB)
                        self.lock.release()

                        lis = out.splitlines()
                        lineacorr = ""
                        for lineacorr in lis:
                            if lineacorr.strip().startswith("Can"):
                                break

                        canI = bool(lineacorr.split()[1])
                        nproc = int(lineacorr.split()[3])

                        # print "I can process:",canI,nproc

                        # nproc = 4 #TEMPORARY FOR TESTING

                        if nproc >= 1:
                            nproc = 1
                            break
                        else:
                            timep.sleep(3)

                while True:
                    # print "ASKED by submit"
                    self.lock.acquire()
                    # print "RECEVEID by submit"
                    out = self.connection.send_command_to_channel(self.channel, cmd, self.promptB)
                    self.lock.release()
                    # scrivere nel Log l'eventuale err.
                    # print "======================OUTPUT===================="
                    # print out
                    # print "======================OUTPUT===================="
                    liout = out.splitlines()
                    entered = False
                    for lineacorr in liout:
                        if lineacorr.strip().startswith("Your"):
                            entered = True
                            break
                    if not entered:
                        print("SGE not available sleeping 10 seconds...")
                        timep.sleep(10)
                        continue
                    nQ = 0
                    #nCluster = int(lineacorr.split()[2])
                    nCluster=lineacorr
                    now = datetime.datetime.now()
                    time = now.strftime("%Y-%m-%d %H:%M")
                    # self.JOBS[job.getName()].append({"timeStart":time,"cmd":cmd,"job":job,"gridCluster":nCluster,"qtotal":1,"queue":nQ})
                    nproc -= 1
                    print("job is: ", nCluster, nQ)
                    break
        except:
            print("Error on submitting a job...")
            print(sys.exc_info())
            traceback.print_exc(file=sys.stdout)
            timep.sleep(3)

    def getStatus(self, jobid, nqueue):
        if jobid not in self.JOBS.keys() and jobid not in self.JOBS_DONE.keys():
            # print "JNP A version"
            return "JNP"  # Job Never Performed
        elif jobid not in self.JOBS.keys() and jobid in self.JOBS_DONE.keys():
            return "NRA"  # Not Registered Anymore
        elif jobid in self.JOBS.keys() and jobid not in self.JOBS_DONE.keys() and len(self.JOBS[jobid]) == 0:
            # print "JNP B version"
            return "JNP"  # Job Never Performed

        while True:
            try:
                cluster = ((self.JOBS[jobid])[nqueue])["gridCluster"]
                nq = ((self.JOBS[jobid])[-1])["queue"]
                break
            except:
                pass

        if nqueue < 0:
            return "NCV"  # Non Corresponding Values
        dizio = self.getGridQueue()
        if cluster not in dizio.keys():
            ntot = ((self.JOBS[jobid])[nqueue])["qtotal"]
            ((self.JOBS[jobid])[nqueue]) = {}
            completed = True
            for tr in self.JOBS[jobid]:
                if len(tr.keys()) > 0:
                    completed = False
                    break
            if completed and len(self.JOBS[jobid]) == ntot:
                self.JOBS_DONE[jobid] = copy.deepcopy(self.JOBS[jobid])
                del self.JOBS[jobid]
            return "NRA"  # Not Registered Anymore
        return ((self.JOBS[jobid])[nqueue], dizio[cluster])

    def getGridQueue(self, cluster="", queue=""):
        if not hasattr(self, "channel"):
            # p = subprocess.Popen(['qstat'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # out, err = p.communicate()
            # out = out.decode("ascii")
            # err = err.decode("ascii")
            # dit = sge.getJobStatsFromUser(getpass.getuser())
            liOut = self.getQueueRawData()
        else:
            self.lock.acquire()
            out = self.connection.send_command_to_channel(self.channel, PATH_REMOTE_SGEPY + ' -qstat ' + self.qname,
                                                          self.promptB)
            self.lock.release()
            liOut = out.splitlines()

        # scrivere nel Log l'eventuale err.
        joc = {}  # (cluster):{owner:user,submitted:time,runtime:time,status:stat,priorized:int,size:kb,cmd:command}
        toStart = False
        for linea in liOut:
            if linea.startswith("-----"):
                toStart = True
                continue
            if toStart:
                dati = (linea.strip()).split()
                if len(dati) < 7:
                    toStart = False
                    continue
                clus = int(dati[0])
                prior = float(dati[1])
                name = dati[2]
                user = dati[3]
                status = dati[4]
                submitted = dati[5] + " " + dati[6]
                # queue = dati[7]
                # master = int(dati[7])
                joc[clus] = {"owner": user, "submitted": submitted, "status": status, "priorized": prior, "name": name}
        return joc

    def isGridAlive(self):
        pass

    def removeJob(self, jobid, nqueue):
        pass

    def removeCluster(self, jobid):
        pass

    def submitJobs(self, job, nqueue):
        if not os.path.exists("./temp/"):
            os.makedirs("./temp/")

        if hasattr(self, "channel"):
            tarro = tarfile.open(os.path.join("./temp/", "transfert.tar.gz"), "w:gz")

            for from_fi, to_fi in self.FILE_TO_COPY:
                tarro.add(from_fi, arcname=to_fi)
            tarro.close()

            buff = self.copy_local_file(os.path.join("./temp/", "transfert.tar.gz"), "transfert.tar.gz")
            buff = self.connection.send_command_to_channel(self.channel, 'tar -zxf ' + str(
                os.path.join(self.get_remote_pwd(), 'transfert.tar.gz')), self.promptB)
            buff = self.remove_remote_file("transfert.tar.gz")
            os.remove(os.path.join("./temp/", "transfert.tar.gz"))
            self.FILE_TO_COPY = []

        while True:
            if len(threading.enumerate()) < 10:
                # print "Starting Thread Job..."
                new_t = SystemUtility.OutputThreading(self.__submitJobs, job, nqueue)
                new_t.start()
                # print "...started!"
                break
            else:
                print("Too many thread", len(threading.enumerate()))
                timep.sleep(3)
        return 0, nqueue

    def __submitJobs(self, job, nqueue):
        workDir = ""
        if hasattr(self, "channel"):
            workDir = self.get_remote_pwd()

        print("WORKING DIRECTORY:", workDir)
        if job.getName() in self.JOBS.keys():
            raise TypeError(
                'Another job is already registered with ' + job.getName() + " and it is impossible to register two jobs with the same name.")

        if not os.path.exists("./grid_jobs/"):
            os.makedirs("./grid_jobs/")

        all_cmds = []
        basecmd = """
#! /bin/bash
# Template for SGE array job for ARCIMBOLDO-BORGES and ARCIMBOLDO
#  Cambiar el valor de wd, el rango es -t 1-n, empieza en 1
# Una vez empezados los parametros de SGE, no poner lineas en blanco o deja de leerlos.
"""

        basecmd += "#$ -S /bin/bash\n"
        if isinstance(job.initdir, str):
            if job.initdir != "":
                basecmd += "#$ -wd " + str(job.initdir) + "\n" + "#$ -N p" + str(
                    job.getName()) + "\n" + "#$ -q " + self.qname + "\n" + "#$ -t 1-" + str(
                    nqueue) + "\n" + "#$ -r y\n#$ -o /dev/null\n#$ -e /dev/null\n\n"
            else:
                basecmd += "#$ -cwd\n" + "#$ -N p" + str(
                    job.getName()) + "\n" + "#$ -q " + self.qname + "\n" + "#$ -t 1-" + str(
                    nqueue) + "\n" + "#$ -r y\n#$ -o /dev/null\n#$ -e /dev/null\n\n"

        elif isinstance(job.initdir, list):
            basecmd += "#$ -N p" + str(
                job.getName()) + "\n" + "#$ -q " + self.qname + "\n" + "#$ -r y\n#$ -o /dev/null\n#$ -e /dev/null\n\n"

        last_dire = ""
        glcount = 1
        for i in range(nqueue):
            cmd = basecmd
            if isinstance(job.initdir, str):
                if i > 0:
                    break
            elif isinstance(job.initdir, list):
                dire = ""
                nume = 0
                summa = 0
                for el in job.initdir:
                    dire, nume = el
                    summa += nume
                    if i < summa:
                        break
                # print "last_dire",last_dire,dire
                if last_dire != dire:
                    cmd += "#$ -wd " + str(dire) + "\n" + "#$ -t " + str(i + 1) + "-" + str(summa) + "\n\n"
                    last_dire = dire
                else:
                    continue

            cmd += "i=$(($SGE_TASK_ID - 1))\n"
            # print "actual cmd"
            # print cmd
            if not hasattr(self, "channel") and job.executable.endswith(".py"):
                cmd += '' + PATH_LOCAL_PYTHON_INTERPRETER + ' ' + job.executable
            elif job.executable.endswith(".py"):
                cmd += '' + PATH_REMOTE_PYTHON_INTERPRETER + ' ' + job.executable
            else:
                cmd += '' + job.executable

            if len(job.args) > 0:
                data = job.getArgs(True)
                data = data.replace("$(Process)", "$i")

                cmd += " " + data

            if len(job.stdIn) > 0:
                data = job.getStdIn(True)
                data = data.replace("$(Process)", "$i")

                cmd += "<" + data

            if len(job.listFilesOut) > 1 and job.getOutput(True, "error") == job.getOutput(True, "output"):
                data = job.getOutput(True, "output")
                data = data.replace("$(Process)", "$i")

                cmd += ">&" + data
            elif len(job.listFilesOut) > 0:
                data = job.getOutput(True, "output")
                data = data.replace("$(Process)", "$i")

                cmd += ">" + data

            if len(job.listFilesOut) > 1 and job.getOutput(True, "error") != job.getOutput(True, "output"):
                data = job.getOutput(True, "error")
                data = data.replace("$(Process)", "$i")

                cmd += "2>" + data

            cmd += '\n'

            f = open("./grid_jobs/" + job.getName() + "_" + str(glcount) + ".sh", "w")
            f.write(cmd)
            f.close()
            all_cmds.append("./grid_jobs/" + job.getName() + "_" + str(glcount) + ".sh")
            glcount += 1

        sent = 0
        tosend = nqueue
        self.JOBS[job.getName()] = []
        nproc = 0
        timep.sleep(60)  # This should avoid NFS delay errors in telemachus
        for cmd in all_cmds:
            try:
                if not hasattr(self, "channel"):
                    if nproc == 0:
                        while True:
                            if float(self.fraction) == 1.0:
                                nproc = tosend
                                break

                            (canI, nproc) = self.canSubmitJobs(getpass.getuser(), qname=self.qname,
                                                              fraction=self.fraction)
                            # nproc = 4 #TEMPORARY FOR TESTING

                            if nproc >= 1:
                                if tosend < nproc:
                                    nproc = tosend
                                break
                            else:
                                timep.sleep(3)
                    while True:
                        p = subprocess.Popen(['qsub', cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        # print cmd.replace('"','').strip()
                        out, err = p.communicate()
                        out = out.decode("ascii")
                        err = err.decode("ascii")
                        # scrivere nel Log l'eventuale err.
                        # print "======================OUTPUT===================="
                        print(out)
                        # print "======================OUTPUT===================="
                        # print "======================ERROR===================="
                        print(err)
                        # print "======================ERROR===================="
                        entered = False
                        liout = out.splitlines()
                        for lineacorr in liout:
                            if lineacorr.strip().startswith("Your"):
                                entered = True
                                break
                        if not entered:
                            print("SGE not available sleeping 10 seconds...")
                            timep.sleep(10)
                            continue
                        nQ = sent
                        #nCluster = int(lineacorr.split()[2])
                        nCluster = lineacorr
                        now = datetime.datetime.now()
                        time = now.strftime("%Y-%m-%d %H:%M")
                        # self.JOBS[job.getName()].append({"timeStart":time,"cmd":cmd,"job":job,"gridCluster":nCluster,"qtotal":nqueue,"queue":nQ})
                        tosend -= 1
                        sent += 1
                        nproc -= 1
                        print("job is: ", nCluster, nQ)
                        break
                else:
                    if nproc == 0:
                        while True:
                            if float(self.fraction) == 1.0:
                                nproc = tosend
                                break

                            self.lock.acquire()
                            out = self.connection.send_command_to_channel(self.channel,
                                                                          PATH_REMOTE_SGEPY + ' -qsub ' + str(
                                                                              self.qname) + " " + str(
                                                                              self.fraction) + "\n", self.promptB)
                            self.lock.release()

                            lis = out.splitlines()
                            lineacorr = ""
                            for lineacorr in lis:
                                if lineacorr.strip().startswith("Can"):
                                    break

                            canI = bool(lineacorr.split()[1])
                            nproc = int(lineacorr.split()[3])

                            # print "I can process:",canI,nproc

                            # nproc = 4 #TEMPORARY FOR TESTING

                            if nproc >= 1:
                                if tosend < nproc:
                                    nproc = tosend
                                break
                            else:
                                timep.sleep(3)
                    # print "Request for lock in submit",cmd
                    print("COPYING", cmd, "INTO", os.path.join(workDir, os.path.basename(cmd)))
                    self.copy_local_file(cmd, os.path.join(workDir, os.path.basename(cmd)), send_now=True,
                                         remote_path_asitis=True)
                    while True:
                        self.lock.acquire()
                        print("EXECUTING", 'nohup qsub ' + os.path.join(workDir,
                                                                        os.path.basename(cmd) + " | tee nohup.out"))
                        out = self.connection.send_command_to_channel(self.channel,
                                                                      'nohup qsub ' + os.path.join(workDir,
                                                                                                   os.path.basename(
                                                                                                       cmd)) + " | tee nohup.out",
                                                                      self.promptB)
                        self.lock.release()
                        # print "Lock free"
                        # scrivere nel Log l'eventuale err.
                        print(cmd)
                        # print "======================OUTPUT===================="
                        print(out)
                        # print "======================OUTPUT===================="
                        entered = False
                        liout = out.splitlines()
                        for lineacorr in liout:
                            if lineacorr.strip().startswith("Your"):
                                entered = True
                                break
                        if not entered:
                            print("SGE not available sleeping 10 seconds...")
                            timep.sleep(10)
                            continue
                        nQ = sent
                        #nCluster = int(lineacorr.split()[2].split(".")[0])
                        nCluster = lineacorr
                        now = datetime.datetime.now()
                        time = now.strftime("%Y-%m-%d %H:%M")
                        # self.JOBS[job.getName()].append({"timeStart":time,"cmd":cmd,"job":job,"gridCluster":nCluster,"qtotal":nqueue, "queue":nQ})
                        tosend -= 1
                        sent += 1
                        nproc -= 1
                        print("job is: ", nCluster, nQ)
                        break
            except:
                print("Error on submitting a job...")
                print(sys.exc_info())
                traceback.print_exc(file=sys.stdout)
                timep.sleep(3)

    def setRequirements(self, requiString):
        pass

    def setMemory(self, memoryString):
        pass

    def setRank(self, rankString):
        pass

    def getCMD(self, jobid):
        try:
            return (self.JOBS[jobid])["cmd"]
        except:
            raise TypeError('The Job is not registered to the GridManager')

    def getQueueRawData(self, qname='all.q'):
        while True:
            try:
                qpipe = subprocess.Popen(['qstat', '-u', '\"', '*', '\"', '-q', qname], stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)
                break
            except:
                timep.sleep(1)
        return qpipe.communicate()[0].splitlines()

    def getQueueSummaryData(self):
        while True:
            try:
                qpipe = subprocess.Popen(['qstat', '-g', 'c', '-ext'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                break
            except:
                timep.sleep(1)
        return qpipe.communicate()[0].splitlines()

    def getUserQueueRawData(self, username, qname=''):
        if len(qname) == 0:
            while True:
                try:
                    qpipe = subprocess.Popen(['qstat', '-u', username], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    break
                except:
                    timep.sleep(1)
        else:
            while True:
                try:
                    qpipe = subprocess.Popen(['qstat', '-u', username, '-q', qname], stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE)
                    break
                except:
                    timep.sleep(1)

        return qpipe.communicate()[0].splitlines()

    def getSGEQueueList(self):
        while True:
            try:
                qpipe = subprocess.Popen(['qconf', '-sql'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                break
            except:
                timep.sleep(1)
        return qpipe.communicate()[0].splitlines()

    def checkQueueName(self, qname):
        if qname not in self.getSGEQueueList():
            print("\nError given queue name does not exist:", qname, "\nThe available queues are:")
            for q in self.getSGEQueueList():
                print("-", q)
            sys.exit('Program terminated due to bad input.\n')

    def getSGEUserData(self, qname='all.q'):
        outlines = self.getQueueRawData(qname)
        userdata = {}
        reg = {}
        reg['run'] = re.compile('R?[rt]')
        reg['held'] = re.compile('^h')
        reg['queued'] = re.compile('qw')
        reg['err'] = re.compile('^E')

        for line in outlines[2:]:
            cols = line.split()
            if len(cols) == 9:
                if cols[3] not in userdata:
                    userdata[cols[3]] = int(cols[8])
                else:
                    userdata[cols[3]] += int(cols[8])
            elif len(cols) == 8 and reg['queued'].match(cols[4]):
                if cols[3] != getpass.getuser():
                    if cols[3] not in userdata:
                        userdata[cols[3]] = int(cols[7])
                    else:
                        userdata[cols[3]] += int(cols[7])

        return userdata

    def getJobStatsFromUser(self, username):
        outlines = self.getQueueRawData()
        reg = {}
        stat = {}

        usreg = re.compile(username)
        reg['run'] = re.compile('R?[rt]')
        reg['held'] = re.compile('^h')
        reg['queued'] = re.compile('qw')
        reg['err'] = re.compile('^E')
        for key in reg:
            stat[key] = 0
        for i in outlines:
            if usreg.search(i):
                cols = i.split()
                for key in reg:
                    if reg[key].match(cols[4]):
                        stat[key] += 1
        return stat

    def getSGEResources(self, qname="all.q"):
        """Return a dictionary containing current SGE resource status (idle, busy)"""
        stats = {'total': 0, 'idle': 0, 'busy': 0, 'error': 0, 'disabled': 0}
        summary = self.getQueueSummaryData()
        for line in summary:
            if qname in line:
                cols = line.split()
                break
        stats['busy'] = int(cols[2])
        stats['idle'] = int(cols[4])
        stats['error'] = int(cols[12])
        stats['disabled'] = int(cols[14])
        stats['total'] = int(cols[5])
        return stats

    def getSGEPoolPercents(self):
        usersdata = self.getSGEUserData()
        sgeslots = self.getSGEResources()['total']
        userpercent = {}
        for user in usersdata:
            userpercent[user] = usersdata[user] * 100 / sgeslots
        return userpercent

    def canSubmitJobs(self, username, qname='all.q', fraction=0.33):
        """ Check SGE queue and determine if user can submit jobs, return a tuple of (True/False, availslots) """
        userdata = self.getSGEUserData(qname)
        sgeslots = self.getSGEResources()
        # print "QUEUE ",qname
        # Free resources need to be computed from the user queue data and the
        # total number of slots, but not the free slots reported by qstat
        idleslots = sgeslots['total'] - sgeslots['error']
        # result = ()
        used = 0
        for user in userdata:
            if user != username:
                # print user," = ",userdata[user]
                idleslots -= userdata[user]
                used += userdata[user]
        # print "available total ",idleslots," total ",sgeslots['total']," used ",used
        if username in userdata:

            if int(idleslots * fraction) > userdata[username] and idleslots > 0:
                result = (True, (int(idleslots * fraction)) - userdata[username])
                return result
            else:
                result = (False, 0)
                return result
        else:
            if idleslots > 0:
                result = (True, (int(idleslots * fraction)))
            else:
                result = (False, 0)
            return result

    def getStats(self, username):
        out = self.getQueueRawData().splitlines()
        count = 0
        running = 0
        waiting = 0
        unknown = 0
        for i in out:
            if username in i:
                count += 1
                if " r " in i:
                    running += 1
                elif " qw " in i:
                    waiting += 1
                else:
                    unknown += 1
        return {"username": username, "jobcount": count, "running": running, "waiting": waiting, "unknown": unknown}

    def getSummary(self):
        out = self.getQueueSummaryData()
        numbers = out[len(out) - 1].split()
        stat = {}
        stat['in use'] = numbers[2]
        stat['available'] = numbers[4]
        stat['on cluster'] = numbers[5]
        stat['error'] = numbers[12]
        stat['disabled'] = numbers[14]
        return stat


class condorManager(Grid):
    """Manager for Condor jobs"""
    universe = ""
    notification = ""
    nice_user = ""
    should_transfer_files = ""
    when_to_transfer_output = ""
    requiString = ""
    memory = ""
    environment = ""
    rank = ""

    def __init__(self, universe="vanilla", notification="error", nice_user="false", should_transfer_files="IF_NEEDED",
                 when_to_transfer_output="ON_EXIT"):
        super(condorManager, self).__init__()
        self.type_grid = "condor"
        self.universe = universe
        self.notification = notification
        self.nice_user = nice_user
        self.should_transfer_files = should_transfer_files
        self.when_to_transfer_output = when_to_transfer_output

    def submitJob(self, job, isthelast=False, forcesubmit=False):
        if not os.path.exists("./temp/"):
            os.makedirs("./temp/")

        if hasattr(self, "channel"):
            if len(self.FILE_TO_COPY) >= self.CUMULATIVE_TRANSFERING or isthelast:
                tarro = tarfile.open(os.path.join("./temp/", "transfert.tar.gz"), "w:gz")

                for from_fi, to_fi in self.FILE_TO_COPY:
                    tarro.add(from_fi, arcname=to_fi)
                tarro.close()

                buff = self.copy_local_file(os.path.join("./temp/", "transfert.tar.gz"), "transfert.tar.gz")
                buff = self.connection.send_command_to_channel(self.channel, 'tar -zxf ' + str(
                    os.path.join(self.get_remote_pwd(), 'transfert.tar.gz')), self.promptB)
                buff = self.remove_remote_file("transfert.tar.gz")
                os.remove(os.path.join("./temp/", "transfert.tar.gz"))
                self.FILE_TO_COPY = []
                for ju in self.JOB_TO_SUBMIT:
                    a = self.submitJob(ju, forcesubmit=True)
                self.JOB_TO_SUBMIT = []
                return 0, 1
            elif not forcesubmit:
                self.JOB_TO_SUBMIT.append(job)
                return 0, 1

        if job.getName() in self.JOBS.keys():
            raise TypeError(
                'Another job is already registered with ' + job.getName() + " and it is impossible to register two jobs with the same name.")

        if not os.path.exists("./grid_jobs/"):
            os.makedirs("./grid_jobs/")

        cmd = ""
        f = open("./grid_jobs/" + job.getName() + ".cmd", "w")
        if not hasattr(self, "channel") and job.executable.endswith(".py"):
            f.write("executable = " + PATH_LOCAL_PYTHON_INTERPRETER + "\n")
            cmd += "executable = " + PATH_LOCAL_PYTHON_INTERPRETER + "\n"
        elif job.executable.endswith(".py"):
            f.write("executable = " + PATH_LOCAL_PYTHON_INTERPRETER + "\n")
            cmd += "executable = " + PATH_REMOTE_PYTHON_INTERPRETER + "\n"
        else:
            f.write("executable = " + job.executable + "\n")
            cmd += "executable = " + job.executable + "\n"

        f.write("universe = " + self.universe + "\n")
        cmd += "universe = " + self.universe + "\n"
        if len(job.stdIn) > 0:
            f.write("input = " + job.getStdIn(False) + "\n")
            cmd += "input = " + job.getStdIn(False) + "\n"
        if self.nice_user != "":
            f.write("nice_user = " + self.nice_user + "\n")
            cmd += "nice_user = " + self.nice_user + "\n"
        if self.notification != "":
            f.write("notification = " + self.notification + "\n")
            cmd += "notification = " + self.notification + "\n"
        # if len(job.listFiles) > 0:
        #       f.write("input = "+job.getInput()+"\n")
        if len(job.listFiles) > 0:
            f.write("transfer_input_files = " + job.getInput(False) + "\n")
            cmd += "transfer_input_files = " + job.getInput(False) + "\n"
        if self.should_transfer_files != "":
            f.write("should_transfer_files = " + self.should_transfer_files + "\n")
            cmd += "should_transfer_files = " + self.should_transfer_files + "\n"
        if self.when_to_transfer_output != "":
            f.write("when_to_transfer_output = " + self.when_to_transfer_output + "\n")
            cmd += "when_to_transfer_output = " + self.when_to_transfer_output + "\n"
        if self.requiString != "":
            f.write("requirements = " + self.requiString + "\n")
            cmd += "requirements = " + self.requiString + "\n"
        if self.memory != "":
            f.write("request_memory = " + self.memory + "\n")
            cmd += "request_memory = " + self.memory + "\n"
        if self.environment != "":
            f.write("Environment = " + self.environment + "\n")
            cmd += "Environment = " + self.environment + "\n"
        if self.rank != "":
            f.write("Rank = " + self.rank + "\n")
            cmd += "Rank = " + self.rank + "\n"
        if job.maxruntime != -1:
            f.write("maxRunTime = " + str(job.maxruntime) + "\n")
            cmd += "maxRunTime = " + str(job.maxruntime) + "\n"
        if job.periodicRemove != "":
            f.write("periodic_remove = " + job.periodicRemove + "\n")
            cmd += "periodic_remove = " + job.periodicRemove + "\n"
        if len(job.args) > 0 or job.executable.endswith(".py"):
            if job.executable.endswith(".py"):
                f.write("Arguments = " + job.executable + " " + job.getArgs(False) + "\n")
                cmd += "Arguments = " + job.executable + " " + job.getArgs(False) + "\n"
            else:
                f.write("Arguments = " + job.getArgs(False) + "\n")
                cmd += "Arguments = " + job.getArgs(False) + "\n"
        if job.initdir != "":
            f.write("initialdir = " + job.initdir + "\n")
            cmd += "initialdir = " + job.initdir + "\n"
        if len(job.listFilesOut) > 0:
            f.write("output = " + job.getOutput(False, "output") + "\n")
            cmd += "output = " + job.getOutput(False, "output") + "\n"
        if len(job.listFilesOut) > 1:
            f.write("error = " + job.getOutput(False, "error") + "\n")
            cmd += "error = " + job.getOutput(False, "error") + "\n"

        f.write("queue " + str(1) + "\n")
        cmd += "queue = " + str(1) + "\n"
        f.close()

        # timep.sleep(60) #This should avoid NFS delay errors in telemachus
        while True:
            try:
                if not hasattr(self, "channel"):
                    p = subprocess.Popen(['condor_submit', "./grid_jobs/" + job.getName() + ".cmd"],
                                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    out, err = p.communicate()
                    out = out.decode("ascii")
                    err = err.decode("ascii")
                    # scrivere nel Log l'eventuale err.
                    # print "======================OUTPUT===================="
                    print(out)
                    # print "======================OUTPUT===================="
                    # print "======================ERROR===================="
                    print(err)
                    # print "======================ERROR===================="
                else:
                    # timep.sleep(5)
                    self.copy_local_file("./grid_jobs/" + job.getName() + ".cmd", job.getName() + ".cmd")
                    out = self.connection.send_command_to_channel(self.channel,
                                                                  'condor_submit ' + os.path.join(self.get_remote_pwd(),
                                                                                                  job.getName() + ".cmd"),
                                                                  self.promptB)
                    # scrivere nel Log l'eventuale err.
                    # print "======================OUTPUT===================="
                    # print out
                    # print "======================OUTPUT===================="
                liOut = out.split("job(s) submitted to cluster")
                liUno = liOut[0].splitlines()
                nQueue = int((liUno[-1]).strip())
                liDue = liOut[1].split(".")
                nCluster = int((liDue[0]).strip())
                now = datetime.datetime.now()
                time = now.strftime("%Y-%m-%d %H:%M")
                self.JOBS[job.getName()] = {"timeStart": time, "cmd": cmd, "job": job, "gridCluster": nCluster,
                                            "queue": nQueue}
                print("job is: ", nCluster, nQueue)
                break
            except:
                print("Error on submitting a job...")
                print(sys.exc_info())
                traceback.print_exc()
                timep.sleep(3)

        return (nCluster, nQueue)

    def getStatus(self, jobid, nqueue):
        if jobid not in self.JOBS.keys() and jobid not in self.JOBS_DONE.keys():
            return "JNP"  # Job Never Performed
        elif jobid not in self.JOBS.keys() and jobid in self.JOBS_DONE.keys():
            nq = (self.JOBS_DONE[jobid])["queue"]
            if nqueue == nq - 1:
                del self.JOBS_DONE[jobid]
            return "NRA"  # Not Registered Anymore

        while True:
            try:
                cluster = (self.JOBS[jobid])["gridCluster"]
                # print "RAM reading: cluster ",cluster
                nq = (self.JOBS[jobid])["queue"]
                break
            except:
                # print "jobid ",jobid,nqueue,"not ready in RAM"
                # print (self.JOBS[jobid])
                pass

        if nqueue < 0 or nqueue >= nq:
            return "NCV"  # Non Corresponding Values
        try:
            dizio = self.getGridQueue(cluster=str(cluster), queue=str(nqueue))
        except:
            return "JNP"
        if (cluster, nqueue) not in dizio:
            if (cluster, nq - 1) not in dizio:
                self.JOBS_DONE[jobid] = copy.deepcopy(self.JOBS[jobid])
                del self.JOBS[jobid]
            return "NRA"  # Not Registered Anymore
        return (self.JOBS[jobid], dizio[(cluster, nqueue)])

    def getGridQueue(self, cluster="", queue=""):
        secarg = ""
        if cluster != "":
            secarg += str(cluster)
            if queue != "":
                secarg += "." + str(queue)

        if not hasattr(self, "channel"):
            p = subprocess.Popen(['condor_q', secarg], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            out = out.decode("ascii")
            err = err.decode("ascii")

        else:
            out = self.connection.send_command_to_channel(self.channel, 'condor_q ' + secarg, self.promptB)

        out = out.decode('unicode')
        # scrivere nel Log l'eventuale err.
        joc = {}  # (cluster,nid):{owner:user,submitted:time,runtime:time,status:stat,priorized:int,size:kb,cmd:command}
        liOut = out.splitlines()
        toStart = False
        for linea in liOut:
            if linea.startswith(" ID"):
                toStart = True
                continue
            if toStart:
                dati = linea.split()
                if len(dati) != 9:
                    toStart = False
                    continue
                clus, nq = dati[0].split(".")
                clus = int(clus)
                nq = int(nq)
                user = dati[1]
                submitted = dati[2] + " " + dati[3]
                runtime = dati[4]
                status = dati[5]
                prio = int(dati[6])
                size = float(dati[7])
                command = dati[8]
                joc[(clus, nq)] = {"owner": user, "submitted": submitted, "runtime": runtime, "status": status,
                                   "priorized": prio, "size": size, "command": command}
        return joc

    def isGridAlive(self):
        pass

    def removeJob(self, jobid, nqueue):
        if jobid not in self.JOBS.keys():
            return "NRA"  # Not Registered Anymore
        cluster = (self.JOBS[jobid])["gridCluster"]
        nq = (self.JOBS[jobid])["queue"]
        if nqueue < 0 or nqueue >= nq:
            return "NCV"  # Non Corresponding Values
        if not hasattr(self, "channel"):
            p = subprocess.Popen(['condor_rm', '' + str(cluster) + '.' + str(nqueue)], stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            out, err = p.communicate()
            out = out.decode("ascii")
            err = err.decode("ascii")
            # scrivere nel Log l'eventuale err.
            # print "======================OUTPUT===================="
            # print out
            # print "======================ERROR===================="
            return err
        else:
            out = self.connection.send_command_to_channel(self.channel,
                                                          'condor_rm ' + '' + str(cluster) + '.' + str(nqueue),
                                                          self.promptB)
            # print "======================OUTPUT===================="
            # print out

    def removeCluster(self, jobid):
        if jobid not in self.JOBS.keys():
            return "NRA"  # Not Registered Anymore
        cluster = (self.JOBS[jobid])["gridCluster"]
        if not hasattr(self, "channel"):
            p = subprocess.Popen(['condor_rm', '' + str(cluster)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            out = out.decode("ascii")
            err = err.decode("ascii")
            # scrivere nel Log l'eventuale err.
            # print "======================OUTPUT===================="
            # print out
            # print "======================ERROR===================="
            return err
        else:
            out = self.connection.send_command_to_channel(self.channel, 'condor_rm ' + '' + str(cluster), self.promptB)
            # print "======================OUTPUT===================="
            # print out

    def submitJobs(self, job, nqueue):
        if not os.path.exists("./temp/"):
            os.makedirs("./temp/")

        if hasattr(self, "channel"):
            tarro = tarfile.open(os.path.join("./temp/", "transfert.tar.gz"), "w:gz")

            for from_fi, to_fi in self.FILE_TO_COPY:
                tarro.add(from_fi, arcname=to_fi)
            tarro.close()

            buff = self.copy_local_file(os.path.join("./temp/", "transfert.tar.gz"), "transfert.tar.gz")
            buff = self.connection.send_command_to_channel(self.channel, 'tar -zxf ' + str(
                os.path.join(self.get_remote_pwd(), 'transfert.tar.gz')), self.promptB)
            buff = self.remove_remote_file("transfert.tar.gz")
            os.remove(os.path.join("./temp/", "transfert.tar.gz"))
            self.FILE_TO_COPY = []

        if job.getName() in self.JOBS.keys():
            self.JOBS[job.getName()] = None
            # raise TypeError('Another job is yet registered with '+job.getName()+" and it is not possible to register two jobs with the same name.")

        if not os.path.exists("./grid_jobs/"):
            os.makedirs("./grid_jobs/")

        cmd = ""
        f = open("./grid_jobs/" + job.getName() + ".cmd", "w")             #NS: write a command file for the condor grid
        f.write("executable = " + job.executable + "\n")
        cmd += "executable = " + job.executable + "\n"
        f.write("universe = " + self.universe + "\n")
        cmd += "universe = " + self.universe + "\n"
        if len(job.stdIn) > 0:
            f.write("input = " + job.getStdIn(True) + "\n")
            cmd += "input = " + job.getStdIn(True) + "\n"
        if self.nice_user != "":
            f.write("nice_user = " + self.nice_user + "\n")
            cmd += "nice_user = " + self.nice_user + "\n"
        if self.notification != "":
            f.write("notification = " + self.notification + "\n")
            cmd += "notification = " + self.notification + "\n"
        # if len(job.listFiles) > 0:
        #       f.write("input = "+job.getInput()+"\n")
        if len(job.listFiles) > 0:                             #job.listFiles contains the extensions to use
            ecco = job.getInput(True)                       #string  "$(Process).hkl, $(Process).pda, $(Process).ent"  ECCO is $(Process).hkl, $(Process).pda, $(Process).ent
            derte = [int(q.split("_")[0]) for q in ecco.split() if "$(Process)" in q and "_fa" not in q and len(q.split("_"))>1]  #, NS: added _fa not in qSHERLOCK, DERTE is []
            if len(derte) > 0:                     #[] in our case so nothing will happen
                valori = list(range(max(derte)+1))
                dizio = {}
                where = {}
                if isinstance(job.initdir, list):
                    for el in job.initdir:
                        for roota, subFoldersa, fileas in os.walk(el[0]):         #NS: look at all files presents in initdir
                            for fileua in fileas:                                     #
                                pdbfa = os.path.join(roota, fileua)
                                if pdbfa.endswith(".pdb") and len(fileua.split("_")) == 2:
                                    a = fileua[:-4].split("_")[1]
                                    b = fileua[:-4].split("_")[0]
                                    if a not in dizio:
                                        dizio[a] = [int(b)]
                                        where[a] = el[0]
                                    else:
                                        dizio[a].append(int(b))
                elif isinstance(job.initdir, str):
                    for roota, subFoldersa, fileas in os.walk(job.initdir):
                        for fileua in fileas:
                            pdbfa = os.path.join(roota, fileua)
                            if pdbfa.endswith(".pdb") and len(fileua.split("_")) == 2:
                                a = fileua[:-4].split("_")[1]
                                b = fileua[:-4].split("_")[0]
                                if a not in dizio:
                                    dizio[a] = [int(b)]
                                    where[a] = job.initdir
                                else:
                                    dizio[a].append(b)

                for key in dizio:
                   # print key,dizio[key]
                    for e in valori:
                        if e not in dizio[key]:
                            try:
                                os.symlink(os.path.join(where[key],"0_"+str(key)+".pdb"), os.path.join(where[key],str(e)+"_"+str(key)+".pdb"))
                            except:
                                pass

            f.write("transfer_input_files = " + job.getInput(True) + "\n")
            cmd += "transfer_input_files = " + job.getInput(True) + "\n"
        if self.should_transfer_files != "":
            f.write("should_transfer_files = " + self.should_transfer_files + "\n")
            cmd += "should_transfer_files = " + self.should_transfer_files + "\n"
        if self.when_to_transfer_output != "":
            f.write("when_to_transfer_output = " + self.when_to_transfer_output + "\n")
            cmd += "when_to_transfer_output = " + self.when_to_transfer_output + "\n"
        if self.requiString != "":
            f.write("requirements = " + self.requiString + "\n")
            cmd += "requirements = " + self.requiString + "\n"
        if self.memory != "":
            f.write("request_memory = " + self.memory + "\n")
            cmd += "request_memory = " + self.memory + "\n"
        if self.environment != "":
            f.write("Environment = " + self.environment + "\n")
            cmd += "Environment = " + self.environment + "\n"
        if self.rank != "":
            f.write("Rank = " + self.rank + "\n")
            cmd += "Rank = " + self.rank + "\n"
        if job.maxruntime != -1:
            f.write("maxRunTime = " + str(job.maxruntime) + "\n")
            cmd += "maxRunTime = " + str(job.maxruntime) + "\n"
        if job.periodicRemove != "":
            f.write("periodic_remove = " + job.periodicRemove + "\n")
            cmd += "periodic_remove = " + job.periodicRemove + "\n"
        if len(job.args) > 0 or job.executable.endswith(".py"):
            if job.executable.endswith(".py"):
                f.write("Arguments = " + job.executable + " " + job.getArgs(True) + "\n")
                cmd += "Arguments = " + job.executable + " " + job.getArgs(True) + "\n"
            else:
                f.write("Arguments = " + job.getArgs(True) + "\n")
                cmd += "Arguments = " + job.getArgs(True) + "\n"
        if len(job.listFilesOut) > 0:
            f.write("output = " + job.getOutput(True, "output") + "\n")
            cmd += "output = " + job.getOutput(True, "output") + "\n"
        if len(job.listFilesOut) > 1:
            f.write("error = " + job.getOutput(True, "error") + "\n")
            cmd += "error = " + job.getOutput(True, "error") + "\n"
        if job.initdir != "" and isinstance(job.initdir, str):
            f.write("initialdir = " + job.initdir + "\n")
            cmd += "initialdir = " + job.initdir + "\n"
        if isinstance(job.initdir, list):
            for el in job.initdir:
                f.write("initialdir = " + el[0] + "\n")
                cmd += "initialdir = " + el[0] + "\n"
                f.write("queue " + str(el[1]) + "\n")
                cmd += "queue = " + str(el[1]) + "\n"
        if isinstance(job.initdir, str):
            f.write("queue " + str(nqueue) + "\n")
            cmd += "queue = " + str(nqueue) + "\n"
        f.close()

        # timep.sleep(60) #This should avoid NFS delay errors in telemachus
        while True:
            try:
                if not hasattr(self, "channel"):
                    p = subprocess.Popen(['condor_submit', "./grid_jobs/" + job.getName() + ".cmd"],
                                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    out, err = p.communicate()
                    out = out.decode("ascii")
                    err = err.decode("ascii")
                    # scrivere nel Log l'eventuale err.
                    # print "======================OUTPUT===================="
                    print(out)
                    # print "======================OUTPUT===================="
                    # print "======================ERROR===================="
                    print(err)
                    # print "======================ERROR===================="
                else:
                    # timep.sleep(5)
                    self.copy_local_file("./grid_jobs/" + job.getName() + ".cmd", job.getName() + ".cmd")
                    out = self.connection.send_command_to_channel(self.channel,
                                                                  'condor_submit ' + os.path.join(self.get_remote_pwd(),
                                                                                                  job.getName() + ".cmd"),
                                                                  self.promptB)
                    # scrivere nel Log l'eventuale err.
                    # print "======================OUTPUT===================="
                    print(out)
                    # print "======================OUTPUT===================="
                liOut = out.split("job(s) submitted to cluster")
                liUno = liOut[0].splitlines()
                nQueue = int((liUno[-1]).strip())
                liDue = liOut[1].split(".")
                nCluster = int((liDue[0]).strip())
                now = datetime.datetime.now()
                time = now.strftime("%Y-%m-%d %H:%M")
                self.JOBS[job.getName()] = {"timeStart": time, "cmd": cmd, "job": job, "gridCluster": nCluster,
                                            "queue": nQueue}
                print("job is: ", nCluster, nQueue)
                break
            except:
                print("Error on submitting a job...")
                # print sys.exc_info()
                # traceback.print_exc()
                timep.sleep(3)

        return (nCluster, nQueue)

    def setRequirements(self, requiString):
        self.requiString = requiString

    def setMemory(self, memoryString):
        self.memory = memoryString

    def setEnvironment(self, environmentString):
        self.environment = environmentString

    def setRank(self, rankString):
        self.rank = rankString

    def getCMD(self, jobid):
        try:
            return (self.JOBS[jobid])["cmd"]
        except:
            raise TypeError('The Job is not registered to the GridManager')


class gridJob:
    """A job for a Grid Manager"""
    nameId = ""
    executable = ""
    initdir = ""
    listFiles = []
    listFilesOut = []
    maxruntime = -1
    periodicRemove = ""
    args = []
    stdIn = []
    ARG_FILE = []
    STDI_FILE = []
    change_FILE = []
    change_FILE_OUT = []

    def __init__(self, nameId):
        # check if nameId is a string, because it must be
        self.nameId = ""
        self.executable = ""
        self.initdir = ""
        self.listFiles = []
        self.listFilesOut = []
        self.maxruntime = -1
        self.periodicRemove = ""
        self.args = []
        self.stdIn = []
        self.ARG_FILE = []
        self.STDI_FILE = []
        self.change_FILE = []
        self.change_FILE_OUT = []
        self.nameId = str(nameId)

    def getName(self):
        return self.nameId

    def setExecutable(self, exe):
        self.executable = str(exe)

    def setInitialDir(self, initdir): #NS: listaDirec is provided as a list of TUPLES (directoy, nb of files) ie: [('/localdata/nicolas/TRIAL_CROSSF/bidon/RUN/8_EXP_LIBRARY/1/./0', 9)]
        if isinstance(initdir, list):
            for s in range(len(initdir)):
                initdir[s] = (str(initdir[s][0]), initdir[s][1])
            self.initdir = initdir         #NS: initdir becomes a list of TUPLE of the form: [('/localdata/nicolas/TRIAL_CROSSF/bidon/RUN/8_EXP_LIBRARY/1/./0', 9)], i.e stays the same
        else:
            self.initdir = str(initdir)

    def addInputFile(self, ifile, couldChange):
        """
        NS comments:
        ifile: extension name, i.e ".hkl"
        couldChange: Bool
        """

        self.listFiles.append(str(ifile))        #NS: listfile was empty so far, append the extension to it
        self.change_FILE.append(couldChange)     #NS: change_FILE was empty so far, append the Boolean to it

        try:
            inde = self.args.index(str(ifile))    #Find the index of the extension if already present in self.arg and add it to self.ARG_FILE
            self.ARG_FILE.append(inde)
        except:
            pass

        try:                                     #Find the index of the extension if already present in self.stdIn and add it to self.ARG_FILE
            inde = self.stdIn.index(str(ifile))
            self.STDI_FILE.append(inde)
        except:
            pass

        if isinstance(self.initdir,str):         #Not the case from A_LITE, self.initdir is a string
            if not os.path.join(self.initdir,ifile):
                open(os.path.join(self.initdir,ifile),"w").close()
            elif isinstance(self.initdir,list) and len(ifile.split("_")) == 2:
                numero = int(ifile.split("_")[1].replace(".pdb",""))
                for z,subdir in enumerate(self.initdir):
                    if numero < (z+1)*subdir[1]: #It might require <= I need to check
                        if not os.path.join(subdir[0],ifile):
                            open(os.path.join(subdir[0],ifile),"w").close()

    def addOutputFile(self, ifile, couldChange):
        self.listFilesOut.append(str(ifile))
        self.change_FILE_OUT.append(couldChange)

        try:
            inde = self.args.index(str(ifile))
            self.ARG_FILE.append(inde)
        except:
            pass

    def setMaxRuntime(self, time):
        self.maxruntime = time

    def setPeriodicRemove(self, time):
        self.periodicRemove = time

    def setOutput(self, outPath):
        self.output = outPath

    def setArguments(self, argsList):
        self.args = argsList                 #sets job.args
        self.ARG_FILE = []
        for inde in range(len(self.args)):   #takes all the arguments form the shelxe line (ex: either file.pda or -m -z etc options, if an argument corresponds to a provided file, put it into  ARG_FILE)
            ar = self.args[inde]
            if ar in self.listFiles:         #NS: adds the index of shelxe argument in self.ARG_FILE
                self.ARG_FILE.append(inde)
            elif ar in self.listFilesOut:
                self.ARG_FILE.append(inde)

    def setStdIn(self, stdIn):
        self.stdIn = stdIn
        self.STDI_FILE = []
        for inde in range(len(self.stdIn)):
            ar = self.stdIn[inde]
            if ar in self.listFiles:
                self.STDI_FILE.append(inde)

    def getInput(self, isInQueue):
        if not isInQueue:
            return ', '.join(self.listFiles)
        else:
            liRe = []
            for tr in range(len(self.listFiles)):         #NS: goes through the list of extensions (.hkl etc), but by INDEX
                fil = self.listFiles[tr]                  #fil = current extension
                change = self.change_FILE[tr]             # True
                if change:
                    lou = fil.split(".")                  # returns for ex ('','hkl'), but will return ('_fa','hkl')
                    if not '_fa' in lou[0]:
                        a = lou[0] + ("$(Process)") + "." + lou[1]   #should be 0.hkl?
                    else:
                        debut=lou[0].replace('_fa','')
                        a = debut+("$(Process)")+ '_fa' + "." + lou[1]
                    liRe.append(a)
                else:
                    liRe.append(fil)
            return ', '.join(liRe)                     #$(Process).hkl, $(Process).pda, $(Process).ent

    def getOutput(self, isInQueue, mode):
        lfo = []
        cfo = []
        if mode == "output":
            lfo = self.listFilesOut[0:1]
            cfo = self.listFilesOut[0:1]
        elif mode == "error":
            lfo = self.listFilesOut[1:2]
            cfo = self.listFilesOut[1:2]

        if not isInQueue:
            return ', '.join(lfo)
        else:
            liRe = []
            for tr in range(len(lfo)):
                fil = lfo[tr]
                change = cfo[tr]
                if change:
                    lou = fil.split(".")
                    a = lou[0] + ("$(Process)") + "." + lou[1]
                    liRe.append(a)
                else:
                    liRe.append(fil)
            return ', '.join(liRe)

    def getArgs(self, isInQueue):
        if not isInQueue:
            return ' '.join(self.args)
        else:
            liRe = []
            for index in range(len(self.args)):
                if index in self.ARG_FILE:
                    fil = self.args[index]
                    indr = self.listFiles.index(fil)
                    change = self.change_FILE[indr]
                    if change:
                        lou = fil.split(".")                  # returns for ex ('','hkl'), but will return ('_fa','hkl')
                        a = lou[0] + ("$(Process)") + "." + lou[1]   #should be 0.hkl?
                        liRe.append(a)
                    else:
                        liRe.append(fil)
                else:                                         #If the argument doesn't correspond to a file
                    if self.args[index] == '_fa':
                        liRe.append(("$(Process)")+"_fa")
                    else:
                        liRe.append(self.args[index])

            return ' '.join(liRe)

    def getStdIn(self, isInQueue):
        if not isInQueue:
            return ' '.join(self.stdIn)
        else:
            liRe = []
            for index in range(len(self.stdIn)):
                if index in self.STDI_FILE:
                    fil = self.stdIn[index]
                    indr = self.listFiles.index(fil)
                    change = self.change_FILE[indr]
                    if change:
                        lou = fil.split(".")
                        a = lou[0] + ("$(Process)") + "." + lou[1]
                        liRe.append(a)
                    else:
                        liRe.append(fil)
                else:
                    liRe.append(self.stdIn[index])
            return ' '.join(liRe)
