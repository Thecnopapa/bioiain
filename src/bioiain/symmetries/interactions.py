import os, json




from .elements import MonomerContact, Monomer
from ..visualisation import pymol


def interactions_per_monomer(crystal, monomer, folder=None):

    print("\n>>", monomer)
    if folder is not None and type(monomer) is str:
        data_path = os.path.join(folder, monomer)

        monomer = Monomer.recover(monomer, data_path=data_path, load_structure=True)

    print(monomer.data["contacts"]["relevant"])

    interactions = monomer.data["contacts"]["relevant"]

    script = pymol.PymolScript(folder=".", name="test", pymol_path="/home/iain/bin/miniconda/bin/pymol")
    script.load(monomer.paths["self"], "monomer")
    contact_folder= monomer.paths["contact_folder"]
    for interaction in interactions:
        interaction_data = json.load(open(os.path.join(contact_folder, interaction+".data.json")))
        print(interaction_data["monomer2"])
        imon = Monomer.recover(data_path=os.path.join(folder, interaction_data["monomer2"]["id"]), load_structure=True)
        script.load(imon.paths["self"], "interacting")
        for c in interaction_data["relevant_contacts"]:
            a1 = c["atom1"]
            a2 = c["atom2"]
            script.line(sele1=f"monomer and c. {a1['chain']} and i. {a1['resn']} and n. CA", sele2=f"interacting and c. {a2['chain']} and i. {a2['resn']} and n. CA")













    script.load(crystal.paths["original"], "original", to="pdb")
    script.cell()
    #script.symmetries()
    #script.group()

    script.write_script()
    script.execute()





    exit()








