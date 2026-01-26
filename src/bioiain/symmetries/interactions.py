import os, json




from .elements import MonomerContact, Monomer







def interactions_per_monomer(monomer, folder=None):

    if folder is not None and type(monomer) is str:
        data_path = os.path.join(folder, monomer)

        monomer = Monomer.recover(monomer, data_path=data_path, load_structure=True)

    print(monomer.data["contacts"]["relevant"])





