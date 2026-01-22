import os
from typing_extensions import Self
from copy import deepcopy

from .operations import *
from ..utilities.logging import log
from ..utilities.maths import *
from ..biopython import Chain




class CrystalElement(Chain):

    def _init(self, *args, **kwargs):
        super()._init(*args, **kwargs)
        self.data["symmetry"] = {
            "operation": None,
            "operation_n": None,
            "is_symmetry": False,
            "positions": None,
            "CoM-frac": None,
            "CoM-orth": None,
            "contacts": None,
        }
        self.data["contacts"] = {
            "all": [],
            "paths": None,
        }
        self.sym_elements = []


    def __repr__(self):
        return "<bi.{} id={} op:{}>".format(
            self.__class__.__name__, self.id, self.data["symmetry"]["operation_n"])

    def __str__(self):
        return "<bi.{} id={} op:{}>".format(
            self.__class__.__name__, self.id, self.data["symmetry"]["operation_n"])



    def generate_symmetries(self, crystal, monomers, ligands, contacts:bool=True,
                            threshold:int|float=15, min_contacts:int=10,
                            contact_method:str="min-contacts") -> list[Self]:
        """
        Generates symmetries of crystal elements from a given symmetry element, and (optionally) calculates contacts
        (theoretically between monomers).
        A distance threshold in Angstroms and a minimum number of contacts is required for contact calculations.
        The symmetry Elements generated contain one structure per operation, which is split across different positions.
        The contacts generated are stored at a monomer level.
        :param crystal: Crystal object.
        :param contacts: Whether to calculate contacts between monomers.
        :param threshold: Threshold in Angstroms to consider a CA-CA contact.
        :param min_contacts: Minimum contacts to consider a monomer-element contact.
        :param contact_method: Method to determine contacts between elements.
        :return: List of generated symmetry Elements
        """

        try:
            log(3, "Generating symmetries ({})".format(self.data["info"]["name"]))
            print(self)
            log(4, self, "CoM:", [round(c) for c in find_com(self.get_atoms())])
            params = crystal.data["params"]
            key = crystal.data["crystal"]["group_key"]
            operations = dictio_space_groups[key]["symops"]

        except KeyError:
            log("warning", f"Symmetry could not be generated for {crystal}")
            return []

        frac_element = entity_to_frac(self, params)
        self.is_frac = True
        frac_element_com = find_com(frac_element.get_atoms())
        self.data["symmetry"]["CoM-frac"] = frac_element_com
        self.data["symmetry"]["CoM-orth"] = find_com(self.get_atoms())
        self.data["symmetry"]["positions"] = []
        if contacts:
            self.data["contacts"]["all"] = []
            self.data["contacts"]["threshold"] = threshold
            self.data["contacts"]["min_contacts"] = min_contacts

        for n, operation in operations.items():

            log(4, "Operation: {}".format(n))
            displaced_element = generate_displaced_copy(frac_element.copy(), key=key, op_n=n).copy()
            displaced_element.data = deepcopy(displaced_element.data)
            displaced_element.data["symmetry"] = {
                "operation_n": n,
                "operation": operation,
                "is_symmetry": True,
                "positions": [],
                "CoM-frac": None,
                "CoM-orth": None,
                "contacts": {}
            }
            displaced_element.data.pop("contacts")

            displaced_element.data["info"]["name"] = "{}_op{}".format(self.data["info"]["name"], n)

            if contacts:
                contact = {m.id:MonomerContact(m, displaced_element,
                                               contact_method=contact_method,
                                               threshold=threshold,
                                               min_contacts=min_contacts) for m in monomers if not (
                    n == 1 and displaced_element.id == m.id
                )}



            for atoms in displaced_element.get_atoms():

                if not atoms.is_disordered() > 0:
                    atoms = [atoms]

                for atom in atoms:

                    deltaX = ((atom.coord[0] - frac_element_com[0]) % 1) - 0.5
                    deltaY = ((atom.coord[1] - frac_element_com[1]) % 1) - 0.5
                    deltaZ = ((atom.coord[2] - frac_element_com[2]) % 1) - 0.5

                    new_coordX = frac_element_com[0] + deltaX
                    new_coordY = frac_element_com[1] + deltaY
                    new_coordZ = frac_element_com[2] + deltaZ
                    new_coord = [new_coordX, new_coordY, new_coordZ]

                    position = [(n_coord - d_coord + 99.5) for n_coord, d_coord in zip(new_coord, atom.coord)]
                    for p in position:
                        assert p % 1 == 0
                    position = tuple([int(p) for p in position])
                    displaced_element.data["symmetry"]["positions"].append(position)

                    if any([p >= 10 for p in position]):
                        print(atom.get_full_id())
                        print("Position:", position)
                        print("Original:", atom.coord)
                        print("Deltas:", deltaX, deltaY, deltaZ)
                        print("New coord:", new_coord)
                        log("error", "Orthogonal position in fractional operation", raise_exception=True)

                    atom.coord = new_coord
                    atom.position = position

                    if contacts and atom.id == "CA":
                        for m in monomers:
                            if n == 1 and displaced_element.id == m.id:
                                continue
                            for a in m.get_atoms():
                                if not a.id == "CA":
                                    continue
                                d = get_fractional_distance(a.coord, atom.coord, self.data["params"])
                                print(atom.parent.id[1], atom.parent.resname, "\t", a.parent.id[1], a.parent.resname, "\t", d, end="\t")

                                if d <= threshold**2:
                                    print("true", end="\r")
                                    contact[m.id].add({
                                        "atom1": a.get_full_id(),
                                        "atom2": atom.get_full_id(),
                                        "distance": d,
                                        "below_threshold": True,
                                        "threshold": threshold,
                                        "position": atom.position,
                                    })
                                else:
                                    print("false", end="\r")

            if contacts:

                for m in monomers:
                    if n == 1 and displaced_element.id == m.id:
                        continue
                    displaced_element.data["symmetry"]["contacts"][m.id] = contact[m.id].data
                    if contact[m.id].check_min_contacts():
                        m.data["contacts"]["all"].append(contact[m.id].data)





            displaced_element.data["symmetry"]["positions"] = list(set(
                displaced_element.data["symmetry"]["positions"]))
            displaced_element.data["symmetry"]["CoM-frac"] = find_com(displaced_element.get_atoms())


            self.data["symmetry"]["positions"].extend(displaced_element.data["symmetry"]["positions"])

            self.sym_elements.append(displaced_element)
            log(5, "Element: {}".format(displaced_element.data["info"]["name"]))
            log(5, "Positions:", displaced_element.data["symmetry"]["positions"])

        if contacts:
            log(4, "Elements in contact: {}".format(len(self.data["contacts"]["all"])))
            for c in self.data["contacts"]["all"]:
                log(5, c["name"])

        self.data["symmetry"]["positions"] = list(set(self.data["symmetry"]["positions"]))
        self.export(structure=False)
        return self.sym_elements







class MonomerContact(object):
    """
    Class to calculate and store element-element contact information.
    :param monomer1: First monomer (certain values, such as the parameter dict will be obtained from this entity).
    :param monomer2: Second monomer.
    :param kwargs: Other parameters to store in the contact such as the threshold used.
    """
    def __init__(self, monomer1, monomer2, **kwargs):
        self.data = {
            "name":"Contact: ({}-{}): Unprocessed".format(monomer1.data["info"]["name"], monomer2.data["info"]["name"]),
            "is_contact": None,
            "positions": [],
            "position": None,
            "kwargs": kwargs,
            "monomer1": {
                "name": monomer1.data["info"]["name"],
                "id": monomer1.id,
                "operation": monomer1.data["symmetry"]["operation_n"],
                "is_symmetry": monomer1.data["symmetry"]["is_symmetry"],
                "position": monomer1.data["symmetry"].get("position", None),
                "CoM-frac": monomer1.data["symmetry"]["CoM-frac"],
            },
            "monomer2": {
                "name": monomer2.data["info"]["name"],
                "id": monomer2.id,
                "operation": monomer2.data["symmetry"]["operation_n"],
                "is_symmetry": monomer2.data["symmetry"]["is_symmetry"],
                "position": monomer2.data["symmetry"].get("position", None),
                "CoM-frac": monomer2.data["symmetry"]["CoM-frac"],

            },
        }
        self.a_a = []


    def add(self, c:dict) -> dict:
        """
        Adds a new atom-atom contact to the list atom-atom contacts. It must be adict with at least the following
        keywords: "is_contact" and "name", additionally a value for "below_threshold" might be required depending on
        the contact method.
        :param c: Dict with the atom-atom contact information.
        """
        self.a_a.append(c)
        return c

    def check_min_contacts(self) -> bool:
        """
        Checks whether this element-element contact meets the required minimum contacts below a certain threshold.
        Requires the min_contacts kwarg on initialisation.
        :return: True if the requirement is met, False otherwise.
        """
        self.data["is_contact"] = False
        pos_dict = {}
        for a in self.a_a:
            if a["below_threshold"]:
                if str(["position"]) in pos_dict:
                    pos_dict[str(["position"])] += 1
                else:
                    pos_dict[str(["position"])] = 1
                if pos_dict[str(["position"])]  >= self.data["kwargs"].get("min_contacts", 1):
                    self.data["is_contact"] = True
                    self.data["positions"].append(a["position"])
                    self.data["name"] = repr(self)
        if self.data["is_contact"]:
            self.data["positions"] = list(set(self.data["positions"]))
            return True
        else:
            self.data["name"] = repr(self)
            return False


    def __repr__(self):
        return "Contact: ({} <-> {}): {}, N:{}, T:{}, min:{}".format(
            self.data["monomer1"]["name"], self.data["monomer2"]["name"],
            self.data["is_contact"], len(self.a_a),
            self.data["kwargs"].get("threshold", None), self.data["kwargs"].get("min_contacts", None))









class Monomer(CrystalElement):
    def _init(self, *args, **kwargs):
        self.paths["export_folder"] += "/monomers"
        self.data["info"]["name"] = self.name().replace("cryst", "mon")
        super()._init(*args, **kwargs)

    def __repr__(self):
        return f"<bi.{self.__class__.__name__} id={self.id} name={self.name()}>"

    def __str__(self):
        return f"<bi.{self.__class__.__name__} id={self.id} name={self.name()}>"



class Ligand(CrystalElement):
    def _init(self, *args, **kwargs):
        self.paths["export_folder"] += "/monomers/ligands"
        self.data["info"]["name"] = self.name().replace("cryst", "lig")
        super()._init(*args, **kwargs)

    def __repr__(self):
        return f"<bi.{self.__class__.__name__} id={self.id} name={self.name()}>"

    def __str__(self):
        return f"<bi.{self.__class__.__name__} id={self.id} name={self.name()}>"














