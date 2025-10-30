import os
from unidecode import unidecode



def clean_string(string, allow=(".", "_"), remove_newlines=True):
    string = unidecode(str(string))
    if remove_newlines:
        string = string.replace("\n", "")
    r = ''.join(e for e in string if e.isalnum() or e in allow)
    return r

def get_digits(string, allow=("."), integer = False):

    try:
        if integer:
            return int(''.join(e for e in unidecode(str(string)) if e.isdigit() or e in allow))
        else:
            return float(''.join(e for e in unidecode(str(string)) if e.isdigit() or e in allow))
    except:
        try:
            if os.environ["BI_VERBOSE"]:
                print("No digits found in: {}".format(string))
                print(''.join(e for e in unidecode(str(string)) if e.isdigit() or e in allow))
        except:
            print("No digits found in: {}".format(string))
            print(''.join(e for e in unidecode(str(string)) if e.isdigit() or e in allow))

        return None