import os, json, math

from ..visualisation import pymol, PymolScript

from .operations import coord_operation_entity, entity_to_frac, entity_to_orth
from ..utilities.sequences import MSA



class InteractionProfile:
    def __init__(self, monomer, save_folder=None, threshold=None, export=True, force=False):
        self.monomer = monomer
        if save_folder is None:
            save_folder = "/tmp/bioiain/interactions"
        self.save_folder = save_folder
        if threshold is None:
            threshold = self.monomer.data["crystal"]["contact_threshold"]
        self.threshold = threshold
        self.labels = None
        self.name = f"{monomer.get_name()}_interactions_T{self.threshold}"




    def generate_labels(self, relative=False, export=True, force=False, dataset=None, msa=None):

        if relative:
            return self._generate_relative_labels(export=True, force=force, dataset=dataset, msa=msa)
        else:
            return self._generate_absolute_labels(export=True, force=force)




    def _generate_relative_labels(self, export=True, force=False, dataset=None, msa=None, similarity=95):
        assert dataset is not None and msa is not None
        similar_ids = msa.get_similar(self.monomer.get_name(), similarity=similarity)
        print(similar_ids)
        exit()





            
    

    def _generate_absolute_labels(self, relative=False, export=True, force=False):
        from .elements import Monomer
        threshold = self.threshold
        #script = PymolScript(self.name, folder=self.save_folder)

        if "interactions" not in self.monomer.data:
            self.monomer.data["interactions"]={
                str(threshold):{
                    "threshold": threshold,
                    "label": None
                }
            }

        if str(threshold) not in self.monomer.data["interactions"]:
            self.monomer.data["interactions"][str(threshold)] = {
                "threshold": threshold,
                "label": None
            }
        if self.monomer.data["interactions"][str(threshold)]["label"] is not None:
            #print("USING PRECALCULATED LABELS")
            return self.monomer.data["interactions"][str(threshold)]["label"]

        #print(f"CALCULATING INTERACTIONS at T= {threshold}")

        interactions = self.monomer.data["contacts"]["relevant"]
        contact_folder = self.monomer.paths["contact_folder"]
        monomer_folder = self.monomer.paths["export_folder"]
        ints = []
        #script.load_entity(self.monomer, "monomer")

        for n, interaction in enumerate(interactions):
            data = json.load(open(os.path.join(contact_folder, interaction + ".data.json")))
            mon1_data = data["monomer1"]
            mon2_data = data["monomer2"]

            reverse = False
            if self.monomer.get_name() == mon2_data["name"]:
                mon1_data, mon2_data = mon2_data, mon1_data
                reverse = True

            mon1 = mon1_data["monomer"]
            mon2 = mon2_data["monomer"]
            operation = mon2_data["operation"]
            position = mon2_data["position"]
            #print("operation:", operation)
            #print(position)

            #print(f"Mon  1: {mon1} / Mon 2: {mon2}")
            imon = Monomer.recover(data_path=os.path.join(monomer_folder, mon2), load_structure=True)
            #print(imon)
            if operation is None:
                #print("ASU", mon1, mon2, imon)
                name = f"interacting_{n}_{imon.id}"
                #script.load_entity(imon, name)

                for c in data["relevant_contacts"]:
                    if c["distance"] > threshold:
                        continue
                    a1 = c["atom1"]
                    a2 = c["atom2"]
                    if reverse:
                        a1, a2 = a2, a1
                    ints.append([a1["resn"], "contact"])
                    group = math.ceil(c["distance"])
                    #script.line(f"int_{group}", sele1=f"monomer and c. {a1['chain']} and i. {a1['resn']} and n. CA",
                    #            sele2=f"interacting_{n} and c. {a2['chain']} and i. {a2['resn']} and n. CA")

            else:
                #print("SYM", mon1, mon2, imon)
                frac = entity_to_frac(imon.copy(), imon.data["params"])
                #for pos in data["positions"]:
                    #pos_str = "_".join([str(p) for p in pos])
                    #name = f"interacting_{n}_{frac.id}_{pos_str}"
                    #if pos is not None:
                        # disp = coord_operation_entity(frac.copy(),
                        #                               key=imon.data["crystal"]["group_key"],
                        #                               op_n=operation,
                        #                               # params=imon.data["params"],
                        #                               distance=pos,
                        #                               )
                        # disp = entity_to_orth(disp, imon.data["params"])
                        # script.load_entity(disp, name)

                for c in data["relevant_contacts"]:
                    if c["distance"] > threshold:
                        continue
                    a1 = c["atom1"]
                    a2 = c["atom2"]

                    if reverse:
                        a1, a2 = a2, a1
                    pos = "_".join([str(p) for p in a2["pos"]])
                    ints.append([a1["resn"], "contact"])
                    #script.color(f"monomer and c. {a1['chain']} and i. {a1['resn']} and n. CA", "red")
                    if mon1 == mon2:
                        f"monomer and c. {a2['chain']} and i. {a2['resn']} and n. CA"
                    group = math.ceil(c["distance"])
                    #script.line(f"int_{group}", sele1=f"monomer and c. {a1['chain']} and i. {a1['resn']} and n. CA",
                    #            sele2=f"interacting_{n}_{frac.id}_{pos} and c. {a2['chain']} and i. {a2['resn']} and n. CA")

                #print()
        atoms = self.monomer.atoms(ca_only=True)
        labels = "N" * len(atoms)
        for i in ints:
            pos, atom = [(n, a) for n, a in enumerate(atoms) if a.resnum == i[0]][0]
            if i[1] == "contact":
                labels = labels[:pos] + "C" + labels[pos + 1:]

        labels = ">" + labels + "<"
        #print(labels)
        self.monomer.data["interactions"][str(threshold)]["label"] = labels
        if export:
            self.monomer.export()
        #script.orient("monomer")
        #script_path = script.write_script()
        #print("pymol script saved at:", script_path)
        self.labels = labels
        return self.labels





class PredictedMonomerContacts(object):
    def __init__(self, monomer, labels, label_to_index):
        self.monomer = monomer
        self.labels = labels
        self.label_to_index = label_to_index

        self._set_bfactors()


    def _set_bfactors(self):
        atoms_by_res = self.monomer.atoms(group_by_residue=True, ca_only=False)
        print(len(self.labels), len(atoms_by_res))
        assert len(self.labels) == len(atoms_by_res)
        for lab, (res, res_atoms) in zip(self.labels, atoms_by_res.items()):

            b = self.label_to_index[lab]
            for atom in res_atoms:
                atom.set_bfactor(b)

    def save_structure(self, folder):
        from ..biopython.imports import write_atoms
        fname = f"{self.monomer.get_name()}_predicted_monomer_contacts"
        path = os.path.join(folder, fname)
        os.makedirs(folder, exist_ok=True)
        path = write_atoms(self.monomer, path)
        return path




