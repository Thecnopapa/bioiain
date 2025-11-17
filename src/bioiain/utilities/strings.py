
from unidecode import unidecode
from ..utilities.logging import log



def clean_string(string:str, allow:list[str]=(".", "_"), remove_newlines:bool=True) -> str:
    """
    Remove unwanted characters from string.
    :param string: String to clean
    :param allow: List of special characters allowed (default ".", "_")
    :param remove_newlines: whether tho remove "\n" (default True)
    :return: Clean string
    """
    string = unidecode(str(string))
    if remove_newlines:
        string = string.replace("\n", "")
    r = ''.join(e for e in string if e.isalnum() or e in allow)
    return r


def get_digits(string:str, allow:list[str]=("."), integer:bool= False):
    """
    Parse digits within a string as int or float.
    :param string: Target string
    :param allow: List of special characters to allow (default ".")
    :param integer: Parse as int, default False -> parsed as float
    :return: Parsed int/float
    """
    try:
        if integer:
            return int(''.join(e for e in unidecode(str(string)) if e.isdigit() or e in allow))
        else:
            return float(''.join(e for e in unidecode(str(string)) if e.isdigit() or e in allow))
    except:
        log("warning", "No digits found in: {}".format(string))
        log(0, ''.join(e for e in unidecode(str(string)) if e.isdigit() or e in allow))
        return None



def string_to_list(line:str, delimiter:str=" ") -> list:
    """
    Split string by the specified delimiter. Empty (or \n) are removed.
    :param line: String to split
    :param delimiter: Delimiter to split at (default " ")
    :return: List of strings
    """
    l = line.split(delimiter)
    nl = []
    for c in l:
        c = c.strip()
        if c == "" or c == "\n":
            continue
        nl.append(c)
    return nl


def add_front_0(string, digits=2, zero = "0"):
    ret = ""
    string = str(string)
    for i in range(digits-len(string)):
        ret += zero
    ret += string
    return ret
