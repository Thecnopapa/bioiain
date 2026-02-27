import os, json, math

from ..visualisation import pymol, PymolScript
from ..visualisation.plots import *

from .operations import coord_operation_entity, entity_to_frac, entity_to_orth
from ..utilities.logging import log
from ..utilities.exceptions import SequenceMissmatchException



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
            labs, n_neighbours = self._generate_relative_labels(export=export, force=force, dataset=dataset, msa=msa, in_lab_var=in_lab_var)
            
            if dual:
                labs = self._generate_dual_labels(labs, export=export, dataset=dataset, force=force)

            return labs, n_neighbours

        else:
            labs = self._generate_absolute_labels(export=export, force=force)
            return labs








    def _generate_dual_labels(self, labs, dataset=None, export=True, force=False, use_classes=True, n_classes=5):


        new_labs = []


        surface_res = self.monomer.get_surface_residues(force=force)

        resnums = self.monomer.atoms(group_by_residue=True).keys()

        #print(type(labs))
        if use_classes:
            discrete_labs = []
            for n, (cont, res) in enumerate(zip(labs, resnums)):
                #print(n, cont, res)
                discrete = int(cont//(1/n_classes))
                #print(discrete)
                discrete_labs.append(discrete)

        #print(labs)
        for res in resnums:
            if int(res) in surface_res:
                if use_classes: new_labs.append("O")
                else: new_labs.append(0)
                #print(res, "outer")
            else:
                if use_classes: new_labs.append("I")
                else:new_labs.append(1)
                #print(res, "inner")


        #print(new_labs)
        try:
            assert len(resnums) == len(new_labs) == len(labs)
        except:
            #print(resnums, "\n", new_labs, "\n", labs)
            print(len(resnums) , len(new_labs) , len(labs))
            raise SequenceMissmatchException()

        if use_classes:
            assert len(new_labs) == len(discrete_labs)
            dual_labels = [f"{nl}-{dl}" for nl, dl in zip(new_labs, discrete_labs)]
            print(json.dumps({k: sum([1 for v in dual_labels if k == v]) for k in sorted(set(dual_labels))}, indent=4))
        else:
            dual_labels = list(zip(labs, new_labs))
        #print(dual_labels)
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

        n_neighbours = len(similar_ids)

        if padding > 0:
            return ([0]*padding + rel_label + [0]*padding, n_neighbours)
        else:
            return (rel_label, n_neighbours) 








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
                log("warning", f"Atom at res: {i[0]} has issues")
                continue
 
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
            if len(label_to_index) < 2:
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



    def _agglomerate(self):
        atomlist = self.monomer.atoms(ca_only=True, group_by_residue=False)
        from sklearn.cluster import AgglomerativeClustering
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        agg = AgglomerativeClustering(n_clusters=None, distance_threshold=25, compute_full_tree=True)
        coords = [a.coord for a in atomlist]
        labels = [a.b for a in atomlist]
        data = [coords, labels]
        agg.fit(coords)
        #print(agg)
        #reduced = agg.transform(coords)
        #[print(x) for x in zip(agg.labels_, coords)]
        #print("...\n")
        #print(reduced)
        #print(reduced.shape)
        print("N CLUSTERS:", agg.n_clusters_)

        labelled_coords = zip(agg.labels_, coords, labels)

        cluster_intensities = {}
        for l, c, b in labelled_coords:
            l = str(int(l))
            if l not in cluster_intensities:
                cluster_intensities[l] = {"i":0, "n":0, "r":0,  "a":0, "center":None,"c":None, "cc":None, "o":None}
            cluster_intensities[l]["i"] += b
            cluster_intensities[l]["n"] += 1



        cmap = mpl.colormaps["plasma"]
        cluster_cmap = mpl.colormaps["gist_rainbow"]
        #print("N COLORS:", cluster_cmap.N)

        for k, v in cluster_intensities.items():
            cluster_intensities[k]["a"] = v["i"] / v["n"]
            cluster_intensities[k]["cc"] = cluster_cmap(round(cluster_cmap.N/agg.n_clusters_*int(k)))

        cluster_intensities = {k:v for k, v in sorted(cluster_intensities.items(), key=lambda x: x[1]["a"], reverse=True)}
        #print(json.dumps(cluster_intensities, indent=4))

        max_i = max(c["a"] for c in cluster_intensities.values())
        for n, (k, v) in enumerate(cluster_intensities.items()):
            cluster_intensities[k]["o"] = n
            cluster_intensities[k]["r"] = v["a"] / max_i
            cluster_intensities[k]["c"] = cmap(round((cluster_intensities[k]["r"]*255/2)+127))

        #print(json.dumps(cluster_intensities, indent=4))


        fig, ax = fig3D()
        for n, (k, v) in enumerate(cluster_intensities.items()):
            cluster_intensities[k]["center"] = np.array([0.,0.,0.])

        for c, b, l in zip(coords, labels, agg.labels_):
            l = str(int(l))
            cluster_intensities[l]["center"] += np.array(c)
            ax.scatter(*c, s=b*100, color=cluster_intensities[l]["c"], edgecolors=cluster_intensities[l]["cc"], alpha=0.5)

        for n, (k, v) in enumerate(cluster_intensities.items()):
            cluster_intensities[k]["center"] /=  cluster_intensities[k]["n"]
            if cluster_intensities[k]["o"] < 5:

                ax.text(*cluster_intensities[k]["center"], v["o"], fontsize=50.**cluster_intensities[k]["r"])

        for n, (k, v) in enumerate(cluster_intensities.items()):
            cluster_intensities[k]["center"] =  list(cluster_intensities[k]["center"])
        #print(json.dumps(cluster_intensities, indent=4))

        fig.show()
        #input("Press Enter to close fig")




    def plot_mpl(self):
        import matplotlib as mpl
        cmap=mpl.colormaps["plasma"]
        print(cmap.N)

        fig, ax = fig3D()

        reslist = self.monomer.atoms(ca_only=True, group_by_residue=True).items()
        pcolor = None
        pcoord=None
        for n, (res, atom_list) in enumerate(reslist):
            for atom in atom_list:
                color_index = round(atom.b*25.5)
                #print(atom.b, color_index, cmap(color_index))
                color = cmap(color_index)
                coord = np.array(atom.coord)
                ax.scatter(*atom.coord, s=100, color=color)
                if n != 0:
                    middle = (coord + pcoord) / 2
                    line1 = list(zip(coord, middle))
                    line2 = list(zip(middle, pcoord))
                    ax.plot(*line1, color=color)
                    ax.plot(*line2, color=pcolor)
                pcolor = color
                pcoord = coord



        fig.show()
        input("Press Enter to close fig")






