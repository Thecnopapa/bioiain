import os, sys, json, xmltodict

import bioiain as bi
import subprocess






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

    def analyse(self, filepath):
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
                raise Exception(f"Interaction XML error: {xml}")
            data["pdb_code"] = pdb_code = xml["pdb_code"].upper()
            data["n_interfaces"] = xml["n_interfaces"]
            xml = xml["interface"]
            for interface in xml:
                interfaces[interface["id"]] = {
                    "info": {k:v for k,v in interface.items() if type(v) not in [list, dict]},
                    "molecules": {}
                    }
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
                raise Exception(f"Assembly XML error: {xml}")
            data["pisa_id"] = xml["name"]
            data["multimeric_state"] = xml["multimeric_state"]
            data["assessment"] = xml["assessment"]
            data["n_assembly_groups"] = xml["total_asm"]
            xml = xml["asm_set"]
            if type(xml) != list:
                xml = [xml]
            for assembly_group in xml:
                assembly_serial = assembly_group["ser_no"]
                #print(">>>> serial group:", assembly_serial)
                #print(assembly_group["assembly"]["serial_no"])
                if type(assembly_group["assembly"]) == dict:
                    asss = [assembly_group["assembly"]]
                elif type(assembly_group["assembly"]) == list:
                    asss = assembly_group["assembly"]
                else:
                    raise Exception(f"Assembly group XML error: {assembly_group['assembly']}")
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

    def run_pisa(self, filepath, interfaces=True, assemblies=True):
        try:
            ccp4_path = os.environ["CCP4"]
        except KeyError:
            bi.log("error", "CCP4 not enabled")
            ccp4_path = None
        if ccp4_path is None:
            return None
        print(" ... CCP4 detected at:", ccp4_path)
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
        except:
            return None

        if interfaces or assemblies:
            int_file, ass_fie = self._pisa_to_xml()
            return self.pisa_id, int_file, ass_fie
        else:
            return self.pisa_id


