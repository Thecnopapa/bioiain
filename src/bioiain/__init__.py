import tempfile, os, sys
sys.path.append(os.path.dirname(__file__))

SUBDIR_NAME="bioiain.d"
TEMP_FOLDER = os.path.join(tempfile.gettempdir(), SUBDIR_NAME)
FD = sys.argv[0]
WD = os.getcwd()

from .utilities import *

log("header", "Initialising Bioiain")
log(1, "SUBDIR:", SUBDIR_NAME)
log(1, "TEMP_FOLDER:", TEMP_FOLDER)
log(1, "FD (file path):", FD)
log(1, "WD (working dir):", WD)
log(1, f"$", *sys.argv)

from .base import *

optional_modules = [".aleph", ".machine", ".tools", ".symmetry"]
for m in optional_modules:
    try:
        __import__(m)
    except exceptions.ModuleLoadError as e:
        log("error", e)
        continue
    except exceptions.ModuleNotEnabled as e:
        log("warning", e)
        continue

del module
del optional_modules





