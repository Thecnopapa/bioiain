import os, sys, json, xmltodict

import subprocess
from ..utilities.logging import log
from ..utilities.exceptions import *


log("header", "Importing PISA module...")
try:
    ccp4_path = os.environ["CCP4"]
except KeyError:
    log("error", "PISA: CCP4 not enabled")
    ccp4_path = None
    raise CCP4NotEnabled("Trying to import the PISA module outside the CCP4 shell")
if ccp4_path is None:
    raise CCP4Error("CCP4 Path not found")
log(1,"CCP4 detected at:", ccp4_path)


class PISA(object):
    def __init__(self,
                 pisa_id="temp",
                 pisa_command="pisa",
                 out_folder="pisa_data",
                 xml_folder="pisa_raw"):
        self.pisa_command = pisa_command
        self.pisa_id = pisa_id
        self.out_folder = out_folder
        self.xml_folder = xml_folder
        self.data = None

    def __getitem__(self, key):
        return self.data.get(key, None)

    def run_pisa(self, filepath, interfaces=True, assemblies=True):

        assert ccp4_path is not None

        filepath = os.path.abspath(filepath)
        cmd = [
            self.pisa_command,
            self.pisa_id,
            "-analyse",
            filepath
        ]
        try:
            print(" ... "+ " ".join(cmd))
            subprocess.run(cmd, check=True)
        except Exception as e:
            raise PISAError(e)

        if interfaces or assemblies:
            int_file, ass_fie = self._pisa_to_xml()
            return self.pisa_id, int_file, ass_fie
        else:
            return self.pisa_id

    def delete(self):
        cmd = [
            self.pisa_command,
            self.pisa_id,
            "-erase"
        ]
        try:
            print(" ... "+ " ".join(cmd))
            subprocess.run(cmd, check=True)
        except Exception as e:
            raise PISAError(e)

    def load(self, data_path):
        with open(data_path, "r") as f:
            self.data = json.load(f)
        return self


    def analyse(self, filepath, force=False):
        if not force:
            fname = os.path.basename(filepath)
            out_path = os.path.join(self.out_folder, f"{fname.split(".")[0]}.pisa.json")
            print(out_path)
            if os.path.exists(out_path):
                self.load(out_path)
                print(" ... reusing pisa output from:", out_path)
                return out_path, self["pdb_code"]
        filepath = os.path.abspath(filepath)
        print(f" ... PISA: analysing: {filepath}")
        self.run_pisa(filepath)
        out, code = self.parse_pisa()
        print(f" ... PISA output for: {code} at: {out}\n")
        return out


    def parse_pisa(self):
        interfaces_path = os.path.join(self.xml_folder, f"{self.pisa_id}.interfaces.xml")
        assemblies_path = os.path.join(self.xml_folder, f"{self.pisa_id}.assemblies.xml")

        # print(" ... " + interfaces_path)
        # print(" ... " + assemblies_path)

        assert os.path.exists(assemblies_path)
        assert os.path.exists(interfaces_path)
        data = {}
        interfaces = {}
        molecules = {}
        with open(interfaces_path) as f:
            print(" ... " + "parsing interfaces")
            xml = xmltodict.parse(f.read())["pdb_entry"]
            if xml["status"] != "Ok":
                raise PISAError(f"Interaction XML error: {xml}")
            data["pdb_code"] = pdb_code = xml["pdb_code"].upper()
            data["n_interfaces"] = int(xml["n_interfaces"])
            if data["n_interfaces"] != 0:
                xml = xml["interface"]
                if type(xml) != list:
                    xml = [xml]
                for interface in xml:
                    try:
                        interfaces[interface["id"]] = {
                            "info": {k:v for k,v in interface.items() if type(v) not in [list, dict]},
                            "molecules": {}
                            }
                    except:
                        print(interface)
                        raise
                    i = interfaces[interface["id"]]
                    for molecule in interface["molecule"]:
                        i["molecules"][molecule["id"]] = {
                            "info": {k:v for k,v in molecule.items() if type(v) not in [list, dict]},
                        }
                        if molecule["chain_id"] not in molecules:
                            molecules[molecule["chain_id"]] = {
                                "id": molecule["id"],
                                "chain_id": molecule["chain_id"],
                                "class": molecule["class"],
                                "residues": {}
                            }
                        m = molecules[molecule["chain_id"]]
                        for n, residue in enumerate(molecule["residues"]["residue"]):
                            if type(residue) == str:
                                continue
                            #print(residue)
                            if residue["ser_no"] not in m["residues"]:
                                m["residues"][residue["ser_no"]] = {
                                    "ser_no": residue["ser_no"],
                                    "name": residue["name"],
                                    "seq_num": residue["seq_num"],
                                    "label_seq_num": residue["label_seq_num"],
                                    "interactions": {}
                                    }
                            solv_en = float(residue["solv_en"])
                            if  solv_en != 0:
                                m["residues"][residue["ser_no"]]["interactions"][interface["id"]] = {
                                    "asa": float(residue["asa"]),
                                    "bsa": float(residue["bsa"]),
                                    "solv_en": solv_en,
                                }
                data["interfaces"] = interfaces
                data["molecules"] = molecules

        assemblies = {}
        with open(assemblies_path) as f:
            print(" ... " + "parsing assemblies")
            xml = xmltodict.parse(f.read())["pisa_results"]
            if xml["status"] != "Ok":
                raise PISAError(f"Assembly XML error: {xml}")
            data["pisa_id"] = xml["name"]
            try:
                data["multimeric_state"] = int(xml["multimeric_state"])
            except:
                raise PISAError(f"multimeric_state not found (probably non-cristallographic)")
            data["assessment"] = xml["assessment"]
            data["n_assembly_groups"] = int(xml["total_asm"])
            is_multimeric = data["multimeric_state"] > 1
            if is_multimeric:
                xml = xml["asm_set"]
            else:
                xml = xml["asu_complex"]
                
            if type(xml) != list:
                xml = [xml]
            for assembly_group in xml:
                if is_multimeric:
                    assembly_serial = assembly_group["ser_no"]
                else:
                    assembly_serial = "0"
                    #print(">>>> serial group:", assembly_serial)
                    #print(assembly_group["assembly"]["serial_no"])
                if type(assembly_group["assembly"]) == dict:
                    asss = [assembly_group["assembly"]]
                elif type(assembly_group["assembly"]) == list:
                    asss = assembly_group["assembly"]
                else:
                    raise PISAError(f"Assembly group XML error: {assembly_group['assembly']}")
                #print_children(assembly_group["assembly"])
                for assembly in asss:
                    #print(">>>> assembly:", assembly["serial_no"])
                    assemblies[assembly["serial_no"]] = {
                        "info": {k:v for k,v in assembly.items() if type(v) not in [list, dict]},
                        "assembly_group": assembly_serial,
                        "interfaces": {}
                    }
                    a = assemblies[assembly["serial_no"]]
                    a["info"]["n_interfaces"] = int(assembly["interfaces"]["n_interfaces"])
                    #print(json.dumps(assemblies[assembly["serial_no"]], indent=4))
                    if a["info"]["n_interfaces"] > 0:
                        #print_children(assembly["interfaces"]["interface"])
                        inters = assembly["interfaces"]["interface"]
                        if type(inters) != list:
                            inters = [inters]
                        for interface in inters:
                            #a["interfaces"][interface["id"]] = {"info": {k:v for k,v in assembly.items() if type(v) not in [list, dict]}}
                            a["interfaces"][interface["id"]] = {"id": interface["id"], "dissociates": interface["dissociates"]}

            data["assemblies"] = assemblies
        os.makedirs(self.out_folder, exist_ok=True)
        out_path = os.path.join(self.out_folder, f"{pdb_code}.pisa.json")
        json.dump(data, open(out_path, "w"), indent=4)
        self.data = data
        return out_path, pdb_code




    def _pisa_to_xml(self, interfaces=True, assemblies=True):
        os.makedirs(self.xml_folder, exist_ok=True)
        cmd = [
            self.pisa_command,
            self.pisa_id,
            "-xml",
        ]
        assemblies_path = None
        interfaces_path = None
        if interfaces:
            interfaces_path = f"{self.xml_folder}/{self.pisa_id}.interfaces.xml"
            cmd_i = cmd + ["interfaces"]
            print(" ... "+ " ".join(cmd_i))
            subprocess.run(cmd_i, stdout=open(interfaces_path, "w"))
        if assemblies:
            assemblies_path = f"{self.xml_folder}/{self.pisa_id}.assemblies.xml"
            cmd_a = cmd + ["assemblies"]
            print(" ... "+ " ".join(cmd_a))
            subprocess.run(cmd_a, stdout=open(assemblies_path, "w"))
        return interfaces_path, assemblies_path



