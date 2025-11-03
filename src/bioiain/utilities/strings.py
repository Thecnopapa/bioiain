
from unidecode import unidecode
from ..utilities.prints import log



def clean_string(string:str, allow:list[str]=(".", "_"), remove_newlines:bool=True):
    string = unidecode(str(string))
    if remove_newlines:
        string = string.replace("\n", "")
    r = ''.join(e for e in string if e.isalnum() or e in allow)
    return r


def get_digits(string:str, allow:list[str]=("."), integer:bool= False):

    try:
        if integer:
            return int(''.join(e for e in unidecode(str(string)) if e.isdigit() or e in allow))
        else:
            return float(''.join(e for e in unidecode(str(string)) if e.isdigit() or e in allow))
    except:
        log("warning", "No digits found in: {}".format(string))
        log(0, ''.join(e for e in unidecode(str(string)) if e.isdigit() or e in allow))

        return None