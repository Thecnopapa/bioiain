import os, json, math

from ..visualisation import pymol, PymolScript

from .operations import coord_operation_entity, entity_to_frac, entity_to_orth




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





    def generate_labels(self, relative=False, export=True, force=False, dataset=None, msa=None, dual=False, in_lab_var="label_path"):

        if relative:
            labs = self._generate_relative_labels(export=export, force=force, dataset=dataset, msa=msa, in_lab_var=in_lab_var)
        else:
            labs = self._generate_absolute_labels(export=export, force=force)

        if dual:
            labs = self._generate_dual_labels(labs, export=export, dataset=dataset, force=force)

        return labs





    def _generate_dual_labels(self, labs, dataset=None, export=True, force=False, use_classes=True, n_classes=5):

        if use_classes:
            assert dataset is not None
            print(dataset)
            print(dataset.data["label_to_index"])
            new_labs = []
            for old_class in dataset.data["label_to_index"]:
                for i in range(n_classes):
                    new_labs.append(f"{old_class}-{i}")
            print(new_labs)
            exit()

        surface_res = self.monomer.get_surface_residues(force=force)

        resnums = self.monomer.atoms(group_by_residue=True).keys()


        for res in resnums:
            if int(res) in surface_res:
                new_labs.append(0)
                #print(res, "inner")
            else:
                new_labs.append(1)
                #print(res, "outer")


        #print(new_labs)
        try:
            assert len(resnums) == len(new_labs) == len(labs)
        except:
            #print(resnums, "\n", new_labs, "\n", labs)
            #print(len(resnums) , len(new_labs) , len(labs))
            raise

        dual_labels = list(zip(labs, new_labs))

        return dual_labels






    def _generate_relative_labels(self, export=True, force=False, dataset=None, msa=None, similarity=95, in_lab_var="label_path"):
        assert dataset is not None and msa is not None
        similar_ids = [self.monomer.get_name()]
        similar_ids.extend(msa.get_similar(self.monomer.get_name(), similarity=similarity))
        #print(similar_ids)
        tmp_fasta = f"/tmp/bioiain/alignments/tmp_alignment.fasta"
        open(tmp_fasta, "w").close()
        os.makedirs(os.path.dirname(tmp_fasta), exist_ok=True)
        simlabels = {}
        for simid in similar_ids:
            lab_path = dataset.embeddings[simid][in_lab_var]
            with open(lab_path, "r") as f:
                simlabels[simid] = f.read().strip()

        padding = dataset.embeddings[simid]["padding"]
        #print(simlabels.keys())
        sim_fasta = msa.msa_fasta
        sim_seqs = sim_fasta.parse(key=simlabels.keys())
        #print(sim_seqs.keys())
        tok_fastas = {}
        for k, v in sim_seqs.items():
            if padding > 0:
                label = simlabels[k][padding:-padding]
            else:
                label = simlabels[k]
            #print(k)
            #print("label:", label)
            #print(v)
            #print("alignment:", v[0])
            replaced = v[0]
            for n, s, in enumerate(v[0]):
                if s == "-":
                    continue
                else:
                    #print(label)
                    replaced = replaced[:n] + label[0] + replaced[n+1:]
                    label = label[1:]
            #print("replaced", replaced, "\n\n")
            assert len(label) == 0
            tok_fastas[k] = replaced

        rel_label = []
        for n, t in enumerate(tok_fastas[self.monomer.get_name()]):
            if t == "-":
                continue
            prob = sum([int(v[n] == "C") for v in tok_fastas.values()]) / len(tok_fastas)
            rel_label.append(prob)
        #print(rel_label)
        #print(len(rel_label))
        if export:
            self.monomer.export()

        if padding > 0:
            return [0]*padding + rel_label + [0]*padding
        else:
            return rel_label








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
            try:
                pos, atom = [(n, a) for n, a in enumerate(atoms) if a.resnum == i[0]][0]
            except:
                print(atoms)
                print(i, i[0])
                [print(n, a, a.resnum, i ) for n, a in enumerate(atoms)][0:20]
                print("...")
                raise

            if i[1] == "contact":
                labels = labels[:pos] + "C" + labels[pos + 1:]

        #labels = ">" + labels + "<"
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
    def __init__(self, monomer, labels, label_to_index=None):
        self.monomer = monomer
        self.labels = labels
        if label_to_index is not None:
            if len(label_to_index) <= 2:
                label_to_index = None
        self.label_to_index = label_to_index

        self._set_bfactors()


    def _set_bfactors(self):
        atoms_by_res = self.monomer.atoms(group_by_residue=True, ca_only=False)
        print(len(self.labels), len(atoms_by_res))
        assert len(self.labels) == len(atoms_by_res)
        for lab, (res, res_atoms) in zip(self.labels, atoms_by_res.items()):

            if self.label_to_index is not None:
                b = self.label_to_index[lab]
            else:
                b = lab
            for atom in res_atoms:
                atom.set_bfactor(b)

    def save_structure(self, folder, extra_name="_predicted_monomer_contacts"):
        from ..biopython.imports import write_atoms
        fname = f"{self.monomer.get_name()}{extra_name}"
        path = os.path.join(folder, fname)
        os.makedirs(folder, exist_ok=True)
        path = write_atoms(self.monomer, path)
        return path




