import tempfile, os

SUBDIR_NAME="bioiain.d"
TEMP_FOLDER = os.path.join(tempfile.gettempdir(), SUBDIR_NAME)
FD = os.path.dirname(__file__)
WD = os.getcwd()
from .utilities.logging import log

log("header", "Initialising Bioiain")
log(1, "SUBDIR:", SUBDIR_NAME)
log(1, "TEMP_FOLDER:", TEMP_FOLDER)
log(1, "FD:", FD)
log(1, "WD:", WD)



__all__ = ["aleph", "base", "machine", "tools", "utilities", "visualisation" ]

__all__.extend(["log", "SUBDIR_NAME", "TEMP_FOLDER", "WD", "FD"])





