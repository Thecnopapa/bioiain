import os, json



from ..visualisation import pymol


from .elements import MonomerContact, Monomer
from .operations import coord_operation_entity, entity_to_frac, entity_to_orth




def get_interaction_profile(monomer, folder):



    interactions = monomer.data["contacts"]["relevant"]
    contact_folder= monomer.paths["contact_folder"]

    for n, interaction in enumerate(interactions):
        data = json.load(open(os.path.join(contact_folder, interaction+".data.json")))
        mon1_data=data["monomer1"]
        mon2_data=data["monomer2"]

        reverse = False
        if monomer.get_name() == mon2_data["name"]:
            mon1_data, mon2_data = mon2_data, mon1_data
            reverse = True

        mon1 = mon1_data["monomer"]
        mon2 = mon2_data["monomer"]
        operation = mon2_data["operation"]
        position = mon2_data["position"]
        print("operation:", operation)
        print(position)

        print(f"Mon  1: {mon1} / Mon 2: {mon2}")
        imon = Monomer.recover(data_path=os.path.join(folder, mon2), load_structure=True)
        print(imon)


        ints= []
        if operation is None:
            print("ASU", mon1, mon2, imon)
            name = f"interacting_{n}"

            for c in data["relevant_contacts"]:
                a1 = c["atom1"]
                a2 = c["atom2"]

                if reverse:
                    a1, a2 = a2, a1
                ints.append([a1["resn"], "contact"])

        else:
            print("SYM", mon1, mon2, imon)
            frac = entity_to_frac(imon.copy(), imon.data["params"])
            for pos in data["positions"]:
                pos_str = "_".join([str(p) for p in pos])
                name = f"interacting_{n}_{pos_str}"
                if pos is not None:
                    disp = coord_operation_entity(frac.copy(),
                                  key=imon.data["crystal"]["group_key"],
                                  op_n=operation,
                                  #params=imon.data["params"],
                                  distance=pos,
                                  )
                    disp = entity_to_orth(disp, imon.data["params"])

                    print(disp, pos)

            for c in data["relevant_contacts"]:
                a1 = c["atom1"]
                a2 = c["atom2"]

                if reverse:
                    a1, a2 = a2, a1
                pos = "_".join([str(p) for p in a2["pos"]])
                ints.append([a1["resn"], "contact"])

            print()
    atoms = monomer.atoms(ca_only=True)
    labels = "N"*len(atoms)
    for i in ints:
        pos, atom = [(n, a) for n, a in enumerate(atoms) if a.resnum == i[0]][0]
        if i[1] == "contact":
            labels = labels[:pos]+"C"+labels[pos+1:]
    print(labels)
    labels = ">"+labels+"<"
    return labels






def interactions_per_monomer(monomer, folder=None, script = None):

    print("\n>>", monomer)
    if folder is not None and type(monomer) is str:
        data_path = os.path.join(folder, monomer)

        monomer = Monomer.recover(monomer, data_path=data_path, load_structure=True)

    print(monomer.data["contacts"]["relevant"])

    interactions = monomer.data["contacts"]["relevant"]
    if script is None:
        script = pymol.PymolScript(name=monomer)
    script.load(monomer.paths["self"], "monomer")
    contact_folder= monomer.paths["contact_folder"]


    embeddings = SaProtEmbeddings(entity=monomer)
    print(embeddings)
    print(embeddings.entity)
    print(embeddings.sequence)

    if embeddings.generate_embeddings() is None: return None
    print(embeddings.embeddings)

    ints = []

    print()
    for n, interaction in enumerate(interactions):
        data = json.load(open(os.path.join(contact_folder, interaction+".data.json")))
        mon1_data=data["monomer1"]
        mon2_data=data["monomer2"]

        reverse = False
        if monomer.get_name() == mon2_data["name"]:
            mon1_data, mon2_data = mon2_data, mon1_data
            reverse = True

        mon1 = mon1_data["monomer"]
        mon2 = mon2_data["monomer"]
        operation = mon2_data["operation"]
        position = mon2_data["position"]
        print("operation:", operation)
        print(position)

        print(f"Mon  1: {mon1} / Mon 2: {mon2}")
        imon = Monomer.recover(data_path=os.path.join(folder, mon2), load_structure=True)
        print(imon)



        if operation is None:
            print("ASU", mon1, mon2, imon)
            name = f"interacting_{n}"
            script.load(imon.paths["self"], name)

            for c in data["relevant_contacts"]:
                a1 = c["atom1"]
                a2 = c["atom2"]

                if reverse:
                    a1, a2 = a2, a1
                ints.append([a1["resn"], "contact"])
                script.line(f"int_{n}", sele1=f"monomer and c. {a1['chain']} and i. {a1['resn']} and n. CA", sele2=f"interacting_{n} and c. {a2['chain']} and i. {a2['resn']} and n. CA")
            script.disable(name)

        else:
            print("SYM", mon1, mon2, imon)
            frac = entity_to_frac(imon.copy(), imon.data["params"])
            for pos in data["positions"]:
                pos_str = "_".join([str(p) for p in pos])
                name = f"interacting_{n}_{pos_str}"
                if pos is not None:
                    disp = coord_operation_entity(frac.copy(),
                                  key=imon.data["crystal"]["group_key"],
                                  op_n=operation,
                                  #params=imon.data["params"],
                                  distance=pos,
                                  )
                    disp = entity_to_orth(disp, imon.data["params"])

                    print(disp, pos)
                    script.load_entity(disp, name)
                    script.disable(name)

            for c in data["relevant_contacts"]:
                a1 = c["atom1"]
                a2 = c["atom2"]

                if reverse:
                    a1, a2 = a2, a1
                pos = "_".join([str(p) for p in a2["pos"]])
                ints.append([a1["resn"], "contact"])
                script.line(f"int_{n}", sele1=f"monomer and c. {a1['chain']} and i. {a1['resn']} and n. CA", sele2=f"interacting_{n}_{pos} and c. {a2['chain']} and i. {a2['resn']} and n. CA")


            print()
    atoms = monomer.atoms(ca_only=True)
    labels = "N"*len(atoms)
    for i in ints:
        pos, atom = [(n, a) for n, a in enumerate(atoms) if a.resnum == i[0]][0]
        if i[1] == "contact":
            labels = labels[:pos]+"C"+labels[pos+1:]
    print(labels)
    labels = "<"+labels+">"
    return labels

    embeddings.add_label_from_string(labels)
    return embeddings

















    exit()








