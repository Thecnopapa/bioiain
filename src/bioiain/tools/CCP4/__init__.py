import os
from ...utilities.exceptions import *

log("header", "Importing CCP4 module...")
try:
    CCP4_PATH = os.environ["CCP4"]
except KeyError:
    log("warning", "CCP4 not enabled")
    CCP4_PATH = None
    raise CCP4NotEnabled("Trying to import the CCP4 module outside the CCP4 shell")
if CCP4_PATH is None:
    raise CCP4Error("CCP4 Path not found")
log(1,"CCP4 detected at:", CCP4_PATH)


__all__ = ["PISA", "TRACER" ,"CCP4_PATH"]
