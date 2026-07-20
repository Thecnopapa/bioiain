import os, sys, shutil, json, requests, time
from io import TextIOWrapper

from . import string_to_list, clean_string
from .logging import log
from .. import WD, SUBDIR_NAME, TEMP_FOLDER
from itertools import accumulate
from ..base import BIEntity

rcsb_pdb_url = "https://files.rcsb.org/download/{}.pdb"
rcsb_cif_url = "https://files.rcsb.org/download/{}.cif"


def relative_path(path, relative_to=None):
    if path is None:
        return None
    if relative_to is None:
        relative_to = WD
    path = os.path.abspath(path)
    relative_to = os.path.abspath(relative_to)
    #print(os.path.commonpath([path, relative_to]))
    if os.path.commonpath([path, relative_to]) in ("", "/"):
        log("warning", "No common path found for {} and {} \nReturning absolute instead...".format(path, relative_to))
    return "./"+os.path.relpath(path, relative_to)




class StructureDataset(object):
    def __init__(self, name="dataset", folder=None, shared_source=True, ignore_blacklist=False):
        self.data = {}
        self.name = name
        self.blacklist = []
        self._blacklist_lock = False
        log(1, "Initialising dataset:", self.name, f"(blacklist={not ignore_blacklist})")

        if folder is None:
            if shared_source:
                folder = os.path.join(SUBDIR_NAME, "data", "shared")
            else:
                folder = os.path.join(SUBDIR_NAME, "data", name)
        self.folder = folder
        os.makedirs(self.folder, exist_ok=True)
        self.blacklist_file = os.path.join(self.folder, f"{self.name}.black.list")
        if os.path.exists(self.blacklist_file) and not ignore_blacklist:
            with open(self.blacklist_file, "r") as bl:
                for line in bl:
                    self.blacklist.append(line.strip().replace("\n", "").split(":")[0].strip())
        else:
            with open(self.blacklist_file, "w") as bl:
                bl.write(f"{self.blacklist_file}\n")
        self.blacklist = list(set(self.blacklist))
        log(2, f"{len(self.blacklist)} paths in blacklist: {self.blacklist_file}")




    class Entry(object):
        def __init__(self, data:dict, dataset):
            self.data = data
            self.dataset = dataset


        def __getattr__(self, key):
            return self.data.get(key)

        def __repr__(self):
            return f"<bi.{self.dataset.__class__.__name__}.{self.__class__.__name__}: {self.code} ({self.name}) at {self.path if self.path is not None else self.url} from {self.source}>"

        def blacklist(self, error=None, reason=None):
            self.dataset.add_to_blacklist(self.path, reason=reason)


    def add_to_blacklist(self, path, error=None, reason=None):
        while self._blacklist_lock:
            print("waiting for blacklist lock...")
            time.sleep(1)
        self._blacklist_lock = True
        try:
            if reason is None:
                reason = error
            reason = str(reason).replace("\n", " /")
            log("warning", f"Blacklisted: {path} ({error.__class__.__name__}):{reason}")
            with open(self.blacklist_file, "a") as bl:
                bl.write(f"{path}: {reason}\n")
            self.blacklist.append(path)
        except:
            self._blacklist_lock = False
            raise
        self._blacklist_lock = False

    def check_blacklist(self, path):
        return path in self.blacklist

    def codes(self) -> list:
        return [e.get("code", None) for e in self.data.values() if not self.check_blacklist(e["path"])]

    def urls(self) -> list:
        return [e.get("url", None) for e in self.data.values() if not self.check_blacklist(e["path"])]

    def paths(self) -> list:
        return [e.get("path", None) for e in self.data.values() if not self.check_blacklist(e["path"])]

    def entities(self, entity_class=BIEntity, return_entries=False, **kwargs):
        for entry in self:
            entity = entity_class.from_file(entry.path, code=entry.name, **kwargs)
            if return_entries:
                yield entity, entry
            else:
                yield entity
            
    def export(self, **kwargs):
        [e.export() for e in self.entities(**kwargs)]

    def load(self, **kwargs):
        [e.atoms() for e in self.entities(**kwargs)]

    def shuffle(self):
        import random

        random_keys = list(self.data.keys())
        random.shuffle(random_keys)
        new_dict = {}
        for k in random_keys:
            new_dict[k] = self.data[k]
        self.data = new_dict



    def __repr__(self):
        return f"<bi.{self.__class__.__name__}: {self.name} N={len(self)}>"

    def __getitem__(self, item):
        return self.Entry(self.data[self.codes()[item]], dataset=self)

    def __len__(self):
        return sum([1 for e in self.data.values() if not self.check_blacklist(e["path"])])

    def __iter__(self):
        self._codes = self.codes()
        self.i = 0
        return self

    def __next__(self):
        if self.i < len(self._codes):
            #c = self._codes[self.i]
            self.i += 1
            try:
                return self[self.i-1]
            except IndexError:
                raise IndexError(self.i -1, len(self._codes))
        else:
            raise StopIteration

    def get(self, code, exception="raise"):
        if exception == "raise":
            return self.data.get(code)
        else:
            return self.data.get(code, exception)


    def add(self, code, name=None, path=None, url=None, source="manual", extension="cif", replace=True, **extras):
        if path in self.blacklist:
            log("warning", f"Path: {path} in blacklist: {self.blacklist_file}")
            return self
        if os.path.getsize(path) > 1.5 * 1024 * 1024:
            log("warning", f"File: {os.path.basename(path)} too large! ({os.path.getsize(path) / 1024 / 1024 :3.2f} MiB)")
            self.add_to_blacklist(path, reason=f"File: {os.path.basename(path)} too large! ({os.path.getsize(path) / 1024 / 1024 :3.2f} MiB)")
            return self
        if code in self.codes():
            if replace:
                log("warning", f"Replacing entry {code} in dataset: {self}")
            else:
                log("warning", f"Ignoring entry {code}, already in dataset: {self}")
                return self

        self.data[code] = {
            "code": code,
            "name": name,
            "path": relative_path(path),
            "url": url,
            "source": source,
            "extension": extension,
            **extras
        }
        return self

    def add_dir(self, folder, exclude:str|list|None=None, replace=True, allowed_extensions=("cif", "pdb"), **extras):
        log(2, f"Adding directory files: {folder} to {self}")
        if exclude is None:
            exclude = []
        elif type(exclude) is str:
            exclude = [exclude]
        counter = 0
        for file in os.listdir(folder):
            if file in exclude:
                continue
            extension = file.split(".")[-1]
            if extension not in allowed_extensions:
                continue
            code = file.split(".")[0]
            name = code
            path = relative_path(os.path.join(folder, file))
            self.add(code, name=name, path=path, url=None, source="dir", extension=extension, replace=replace, **extras)
            counter += 1
        log(2, f"Added {counter} files to {self}")

        return self

    @classmethod
    def from_dir(cls, folder, name=None, exclude:str|list|None=None, replace=True, ignore_blacklist=False, **extras):
        if name is None:
            name = os.path.dirname(os.path.abspath(folder))
        self = cls(name=name, ignore_blacklist=ignore_blacklist)
        self.add_dir(folder, exclude=exclude, replace=replace,**extras)
        return self


    def add_list(self, file_or_list:str|list, download_as="cif", base_url=None, force=False, replace=True, **extras):
        log(2, f"Adding list: {file_or_list} to {self}")
        pdb_links = []
        pdbs_to_download = []
        pdb_paths = []
        if type(file_or_list) is list:
            f = sorted(file_or_list)
        else:
            f = open(file_or_list)

        for line in f:
            line = line.split("#")[0].replace("\n", "")
            new = string_to_list(line, delimiter=",")

            for n in new:
                if n.strip == "":
                    continue
                if "http" in n:
                    pdb_links.append(n)
                elif os.path.exists(n):
                    pdb_paths.append(n)
                else:
                    n = clean_string(n.split(".")[0])
                    pdbs_to_download.append(n.upper())
        if type(f) is TextIOWrapper:
            f.close()

        if base_url is None:
            if download_as.lower() == "pdb":
                base_url = rcsb_pdb_url
            elif download_as.lower() == "cif":
                base_url = rcsb_cif_url
        url_file = os.path.join(TEMP_FOLDER, "urls")
        os.makedirs(url_file, exist_ok=True)
        url_file = os.path.join(url_file, f"{self.name}.url.list")
        with open(url_file, "w") as f:
            for url in sorted(pdb_links):
                f.write(url+"\n")
            for code in sorted(pdbs_to_download):
                url = base_url.format(code)
                f.write(url+"\n")

        with open(url_file) as lf:
            counter = 0
            failed_counter = 0
            skipped_counter = 0
            for line in lf:
                line = line.replace("\n", "")
                f_name = line.split("/")[-1]
                subfolder = os.path.join(self.folder, f_name[0])
                os.makedirs(subfolder, exist_ok=True)
                f_path = os.path.join(subfolder, f_name)
                code = f_name.split(".")[0]
                extension = f_name.split(".")[-1]
                url = line
                if not os.path.exists(f_path) or force:
                    log(3, f"Downloading {url}...", end="\r")
                    response = requests.get(url)
                    if response.status_code != 200:
                        log("Error", "Failed to download from:", line)
                        failed_counter += 1
                    else:
                        with open(f_path, "w") as f:
                            f.write(response.text)
                        counter += 1
                else:
                    skipped_counter += 1
                self.add(code, name=code, path=f_path, url=url, source="list_downloaded", extension=extension, replace=replace,**extras)
        if counter > 0 or failed_counter > 0:
            print()
        log(2, f"{counter} files downloaded, {failed_counter} failed, {skipped_counter} skipped")
        for path in pdb_paths:
            code = os.path.basename(path.split(".")[0])
            extension = f_name.split(".")[-1]
            self.add(code, name=code, path=path, url=None, source="list_path", extension=extension, replace=replace, **extras)
        return self

    @classmethod
    def from_list(cls, file_or_list:str|list, name=None, download_as="cif", base_url=None, force=False, replace=True, ignore_blacklist=False, **extras):
        if name is None:
            if type(file_or_list) is list:
                name = "pdb_list"
            elif type(file_or_list) is str:
                name = os.path.basename(file_or_list).split(".")[0]
        self = cls(name=name, ignore_blacklist=ignore_blacklist)
        self.add_list(file_or_list, download_as=download_as, base_url=base_url, force=force, replace=replace, **extras)
        return self
