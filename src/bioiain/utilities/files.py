import os, sys, shutil, json
from .logging import log




def relative_path(path, relative_to=None):
    if path is None:
        return None
    if relative_to is None:
        relative_to = os.abspath(os.getcwd())
    path = os.path.abspath(path)
    relative_to = os.path.abspath(relative_to)
    #print(os.path.commonpath([path, relative_to]))
    if os.path.commonpath([path, relative_to]) in ("", "/"):
        log("warning", "No common path found for {} and {} \nReturning absolute instead...".format(path, relative_to))
    return "./"+os.path.relpath(path, os.path.dirname(relative_to))







