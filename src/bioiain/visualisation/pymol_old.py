import os

from ..utilities import *
from ..utilities.maths_old import *

import pymol


colours = ['green', 'cyan', 'red', 'yellow', 'violet','blue',
               'salmon', 'lime', 'pink', 'slate', 'magenta', 'orange', 'marine',
               'olive', 'purple', 'teal', 'forest', 'firebrick', 'chocolate',
               'wheat', 'white', 'grey']
ncolours = len(colours)


mpl_colours = 'blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan'
mpl_ncolours = len(mpl_colours)

def pymol_colour_everything(start_at=0):
    # List of available colours

    num_states = pymol.cmd.count_states("(all)")
    # Loop over objects
    i = 0
    for obj in get_all_obj():
        for state in range(1,num_states):
            pymol.cmd.color(colours[i+start_at], obj+" and state "+str(state))
            i += 1
            if i >= ncolours:
                i = 0

def pymol_start(show=False, quiet=True):
    # Start PyMol

    print("(PyMol) Starting PyMol...")
    if show:  # Debug stuff / -c: no interface, -q no startup message, -Q completely quiet
        pymol.finish_launching(["pymol"])  #
        if quiet:
            pymol.finish_launching(["pymol", "-q"])
    else:
        pymol.finish_launching(["pymol", "-cqQ"])


def pymol_close():
    print("(PyMol) Closing PyMol...")
    pymol.cmd.quit()

def pymol_reset():
    pymol.cmd.delete("All")

def pymol_reinitialize():
    pymol.cmd.reinitialize()



def pymol_delete(sele = "all", quiet=False):
    if not quiet:
        print("(PyMol) Deleting:", sele)
    pymol.cmd.delete(name = sele)


def pymol_hide(sele, hide = None, silent = True):
    all_reprs = ["lines","spheres","mesh","ribbon","cartoon","sticks","dots","surface","labels","nonbonded","nb_spheres"]
    if hide is not None and len(sele) >0:
        if hide == "all":
            hide = " ".join(all_reprs)
        if not silent:
            print("(PyMol) Hiding:", sele)
        pymol.cmd.hide(representation=hide, selection=sele)

def pymol_format(representation,identifier="", hide="all", colour = None, spectrum=None, quiet=False):
    for obj in get_all_obj():
        if identifier in obj:
            pymol_hide(obj, hide=hide)
            if not quiet:
                print("(PyMol) Formatting:", obj, representation)
            pymol.cmd.show(representation=representation, selection=obj)
            if colour is not None:
                pymol_colour(colour, obj=obj, spectrum=spectrum)


def pymol_colour(colour, obj = "(all)", sele = None, spectrum=None, silent =False, minimum = None, maximum=None):
    if sele is not None:
        sele_str = "({} and {})".format(obj,sele)
    else:
        sele_str = "({})".format(obj)
    if not silent:
        print("(PyMol) Colouring:", sele_str, colour, spectrum)
    if spectrum is not None:
        pymol.cmd.spectrum(spectrum, colour, sele_str, minimum=minimum, maximum=maximum)
    elif colour == "chainbow":
        palette = "rainbow"
        if spectrum is not None:
            palette = spectrum
        pymol.cmd.util.chainbow(sele_str, palette=palette)
    else:
        pymol.cmd.colour(colour, sele_str)


def pymol_save_small(file_name, folder, dpi=300, height=100, width=150, save_session=None):
    print("(PyMol) Saving small:", folder, file_name)
    image_path = os.path.join(folder, file_name+".png")
    pymol.cmd.png(image_path, width=width, height=height, dpi=dpi)
    return image_path

def pymol_save_snapshot(file_name, folder, dpi=300, height=400, width=600, orient=True):
    print("(PyMol) Saving snapshot:", folder, file_name)
    image_path = os.path.join(folder, file_name)
    if image_path[-1].isdigit():
        image_path += "ss.png"
    else:
        image_path += ".png"
    #print(os.listdir(folder))
    if orient:
        pymol_orient()
    print(image_path)
    try:
        pymol.cmd.png(image_path, width=width, height=height, dpi=dpi, quiet=0, ray=0)
    except:
        print("(PyMol) Failed to save image:", image_path)
        image_path = None
    #pymol.cmd.png("test.png")
    #pymol.cmd.save(image_path, quiet=0)
    return image_path

def pymol_save_session(file_name, folder, mode="pse"):
    path = os.path.join(folder, file_name + "."+mode)
    pymol.cmd.save(path)
    return path

def get_all_obj():
    return pymol.cmd.get_names(type='objects')



def pymol_align_all(all_obj = None):
    if all_obj is None:
        all_obj = pymol.cmd.get_names(type='objects')
    #print(all_obj)
    obj1 = all_obj[0]
    for obj2 in all_obj:
        if obj2 != obj1:
            pymol_align__obj(obj1, obj2)
    pymol.cmd.orient("(all)")


def pymol_align__obj(obj1, obj2, orient = True, quiet=False, target_state=0, mobile_state=0):
    if quiet:
        quiet = 1
    else:
        quiet = 0
    r = pymol.cmd.align(obj2, obj1, target_state=target_state, mobile_state=mobile_state, quiet=quiet)
    if not quiet:
        print3("aligned:", obj1, obj2)
    if orient:
        pymol.cmd.orient("(all)")
    return r




def pymol_align_chains(chains_to_align, target_state= 0, mobile_state=0):
    print("(PyMol) Aligning chains")
    #all_obj = pymol.cmd.get_names(type='objects')
    #print1(all_obj)
    obj1, chain1 = chains_to_align[0]
    sele1 = "({} and c. {})".format(obj1, chain1)
    for obj2, chain2 in chains_to_align[1:]:
        if obj2 != obj1:
            sele2 = "({} and c. {})".format(obj2, chain2)
            #print(sele1)
            #print(sele2)
            i = 0
            while i <= 5:
                if i >= 5:
                    print(pymol.cmd.get_names(type='objects'))
                    print(sele1, "//", obj1, "//", chain1)
                    print(sele2, "//", obj2, "//", chain2)
                    print("Aligning failed with 5 tries")
                    i += 1
                    continue
                try:
                    pymol_align__obj(sele1, sele2, target_state=target_state, mobile_state=mobile_state)
                    i += 1
                    break
                except:
                    print("Aligning..., try n={}".format(i))
                    print(pymol.cmd.get_names(type='objects'))
                    print(sele1, "//", obj1, "//", chain1)
                    print(sele2, "//", obj2, "//", chain2)
                    session_path = pymol_save_temp_session()
                    pymol_open_session_terminal(session_path)
                    pymol_align__obj(sele1, sele2, target_state=target_state, mobile_state=mobile_state)
                    quit()
                    i += 1



def pymol_align_chains_best(chains_to_align, double_best = False, cluster = None):
    #all_obj = pymol.cmd.get_names(type='objects')
    #print1(all_obj)
    obj1, chain1a, chain1b = chains_to_align[0]
    sele1a = "{} and c. {}".format(obj1, chain1a)
    sele1b = "{} and c. {}".format(obj1, chain1b)

    progress = ProgressBar(len(chains_to_align), silent=True)
    for obj2, chain2a, chain2b in chains_to_align[1:]:
        sele2a = "{} and c. {}".format(obj2, chain2a)
        sele2b = "{} and c. {}".format(obj2, chain2b)
        if obj2 != obj1:
            ali2a = pymol_align__obj(sele1a, sele2a, orient=False, quiet=True)
            ali2a = pymol_align__obj(obj1, obj2, orient=False, quiet=True)
            #print4(ali2a)
            #input()
            ali2b = pymol_align__obj(sele1a, sele2b, orient=False, quiet=True)
            #ali2b = pymol_align__obj(obj1, obj2, orient=False, quiet=True)
            #print4(ali2b)
            #input()
            if ali2a[0] > ali2b[0]:
                ali1a = pymol_align__obj(sele1a, sele2a, orient=False, quiet=True)
                ali1a = pymol_align__obj(obj1, obj2, orient=False, quiet=True), sele2a
            else:
                ali1a = ali2b, sele2b

            if double_best:
                ali2a = pymol_align__obj(sele1b, sele2a, orient=False, quiet=True)
                ali2a = pymol_align__obj(obj1, obj2, orient=False, quiet=True)
                #print4(ali2a)
                #input()
                ali2b = pymol_align__obj(sele1b, sele2b, orient=False, quiet=True)
                ali2b = pymol_align__obj(obj1, obj2, orient=False, quiet=True)
                #print4(ali2b)
                #input()
                if ali2a[0] < ali2b[0]:
                    ali1b = pymol_align__obj(sele1b, sele2a, orient=False, quiet=True)
                    ali1b = pymol_align__obj(obj1, obj2, orient=False, quiet=True), sele2a
                else:
                    ali1b = ali2a, sele2b

                if ali1a[0][0] < ali1b[0][0]:
                    ali = pymol_align__obj(sele1a, ali1a[1], orient=False,quiet=True)

            ali = pymol_align__obj(obj1, obj2, orient=False,quiet=True)
        progress.add(info=cluster)








def pymol_set_state(state):
    print("(PyMol) Setting state to: ", state, )
    pymol.cmd.set("state", state)


def pymol_orient(enabled = True):

    if enabled:
        print("(PyMol) Orienting enabled")
        pymol.cmd.orient("(enabled)")
    else:
        print("(PyMol) Orienting all")
        pymol.cmd.orient("(all)")





def pymol_draw_line(coord1, coord2, name = "d", state = -1, quiet= True):
    if not quiet:
        print("(PyMol) ({}) Line between:".format(name), coord1, "and", coord2, end="\r")
    pymol.cmd.pseudoatom("tmp1", pos=coord1, state=state)
    pymol.cmd.pseudoatom("tmp2", pos=coord2, state=state)
    n=0
    #print("Waiting for pymol... {} currently: {}".format(n, pymol_get_all_objects()), end="\r")

    while "tmp1" in pymol_get_all_objects():
        n+=1
        try:
            pymol.cmd.distance(name, "(tmp1)","tmp2", state=state)
            pymol.cmd.delete("tmp1")
            pymol.cmd.delete("tmp2")
        except:
            pymol.cmd.delete("tmp1")
            pymol.cmd.delete("tmp2")
            pymol.cmd.pseudoatom("tmp1", pos=coord1, state=state)
            pymol.cmd.pseudoatom("tmp2", pos=coord2, state=state)
        print("Waiting for pymol... {}".format(n), end="\r")
def pymol_paint_contacts(obj, contact_list, colour ="yellow", quiet = True):
    print("(PyMol) Colouring contacts in {}".format(obj))
    for chain, resn, in contact_list:
        sele = "c. {} and i. {}".format(chain, resn)
        pymol_colour(colour, obj, sele, silent = quiet)


def pymol_temp_show(structure, disable = False, delete=False, name=None, folder = None):
    from Bio.PDB import PDBIO

    exporting = PDBIO()
    exporting.set_structure(structure)
    all_obj = [n.upper() for n in pymol_get_all_objects()]
    #print(all_obj)
    if name is None:
        name = "temp_{}.pdb".format(structure.id)
        if name.upper() in all_obj:
            n = 1
            while name.upper() in all_obj:
                name = "temp_{}_{}.pdb".format(structure.id, n)
                n += 1
    else:
        if not name.endswith(".pdb"):
            name += ".pdb"
    if folder is None:
        folder = "./temp"
    path = os.path.join(folder, name)
    exporting.save(path)
    pymol_start(show=True)
    if disable:
        pymol_disable("(all)")
    loaded = pymol_load_path(path)
    if delete:
        os.remove(path)
    return loaded

def pymol_get_all_objects():
    return pymol.cmd.get_names(type='objects')


def pymol_move(sele, distance = [50, 0, 0], camera = 0, state = 0):
    pymol.cmd.translate(distance, sele , state,camera)


def pymol_sphere(coords, name = None, colour="white", state = -1, scale = 8):
    if name is None:
        n_spheres = sum(1 for a in ("tmp_sphere" in obj for obj in pymol_get_all_objects()) if a)
        name = "tmp_sphere_{}".format(n_spheres)
    #print(name)
    #pymol.cmd.pseudoatom(pos = coords, object = name)
    pymol.cmd.pseudoatom(object=name, pos=coords, color = colour, elem="Ca", state = state)
    pymol.cmd.set("sphere_scale", scale)#, selection="({})".format(name))
    pymol.cmd.set("sphere_transparency", 0.25)
    pymol_format("sphere", name, quiet= True)







def pymol_paint_single_face(obj_list, face1, face2, color_dict):
    if color_dict is not None:
        color1 = color_dict[face1]
        color2 = color_dict[face2]
    else:
        color1 = mpl_colours[int(face1)%mpl_ncolours]
        color2 = mpl_colours[int(face2)%mpl_ncolours]
    print(obj_list)
    for n, (obj, o) in enumerate(zip(pymol_get_all_objects(), obj_list.values())):
        print(o)
        chain1 = o[1]
        chain1ischain1 = True


        sele1 = "c. {}".format(chain1)
        sele2 = "not c. {}".format(chain1)
        if not chain1ischain1:
            sele1, sele2 = sele2, sele1
        pymol_colour(color1, obj, sele1, silent = False)
        pymol_colour(color2, obj, sele2, silent = False)



def pymol_paint_all_faces(obj, face_dict=None):
    print("(PyMOL) paining all faces in {}".format(obj))
    if face_dict is None:
        from faces import GR_colours, GR_dict
        assert obj.best_fit == "GR"
        for o in pymol_get_all_objects():
            if obj.id in o:
                for face, ress in GR_dict.items():
                    sele = "({})".format(" or ".join("i. {}".format(res) for res in ress))
                    # print(sele)
                    pymol_colour(GR_colours[face], o, sele, silent=True)
    else:
        for o in pymol_get_all_objects():
            if obj in o:
                for face, ress in face_dict.items():
                    sele = "({})".format(" or ".join( "i. {}".format(res) for res in ress))
                    #print(sele)
                    pymol_colour(mpl_colours[int(face)%mpl_ncolours], o, sele, silent=True)



def pymol_save_temp_session(path=None, name="temp_session.pse"):
    if path is None:
        path = os.path.join("./temp", name)
    pymol.cmd.save(path)
    return path

def pymol_open_session_terminal(path):
    open_file_from_system(path)


def pymol_command_in_new_process(path, command):
    import subprocess
    print("Opening:", path)
    subprocess.Popen(["nohup", "xdg-open", path, "&", command], start_new_session=True)

def pymol_save_cluster(obj_list, name="CLUSTER_X.pdb", folder=None, state=0):
    if folder is None:
        folder = "./temp"
    objects = []
    for obj in obj_list:
        if type(obj) is list or type(obj) is tuple:
            obj = obj[0]
        objects.append(obj)
        #print(obj)
    #print(objects)
    path = os.path.join(folder, name)
    pymol.cmd.multisave(path, pattern=" or ".join(objects), state=state)
    return path

def pymol_open_saved_cluster(path, name_list=None, only_even=True, spheres = False):
    from Bio.PDB import PDBParser
    from faces import find_com, get_pca

    cluster_data = {"names": name_list,
                    "models": [],
                    "chains": [],
                    "pcas": [],
                    "path": path,
                    "corners": [],
                    "coms": []}

    cluster_structure = PDBParser(QUIET=True).get_structure(os.path.basename(path.split(".")[0]), path)
    models = cluster_structure.get_models()
    if only_even:
        models = [model for model in models if model.id % 2 == 0]
    spheres = []
    for model, name in zip(models, cluster_data["names"]):
        print(model)
        cluster_data["models"].append(model)
        for chain in model.get_chains():
            print2(chain)
            cluster_data["chains"].append(chain)
            com = find_com(chain.get_atoms())

            pca = get_pca(chain, com=com, closer_to="N")
            cluster_data["coms"].append(com)
            cluster_data["pcas"].append(pca)
            if spheres:
                pymol_sphere(com, name=name + "_" + chain.id + "_com")
            components = pca.components_
            variances = pca.explained_variance_
            ratios = pca.explained_variance_ratio_
            if "inverse" in pca.__dict__.keys():
                if pca.inverse:
                    components[1], components[2] = components[2].copy(), components[1].copy()
                    variances[1], variances[2] = variances[2].copy(), variances[1].copy()
            sphere_coords = com
            corner = [0, 0, 0]
            for n, (component, variance, ratio) in enumerate(zip(components, variances, ratios)):
                #print("sphere-coords:", sphere_coords)
                #print(com, component, variance)
                pymol_draw_line(com, tuple([c + (co * variance) for c, co in zip(com, component)]),
                                name="{}_{}_component_{}".format(name, chain.id, n), quiet=True)
                sphere_coords = add(sphere_coords, tuple([co * variance for co in component]))
                corner = add(corner, tuple([co * variance for co in component]))
            #print("sphere-coords:", sphere_coords)
            spheres.append(sphere_coords)
            cluster_data["corners"].append(corner)
            if spheres:
                pymol_sphere(sphere_coords, name=name+"_"+chain.id + "_corner")
        pymol_group(name, name="-"+name)

def pymol_list_to_bfactors(val_list, obj_name, resids, paint=True):
    print("(PyMol) Altering B-factors:", obj_name)
    for n, val in zip(resids,val_list):
        sele = "{} and i. {}".format(obj_name, n)
        if val == 666:
            val = -1
        #print(val, "-->", sele)
        pymol.cmd.alter(selection=sele, expression="b={}".format(val))

def pymol_split_states(obj_name, start=1, end = 0):
    pymol.cmd.split_states(obj_name, first=start, last=end)
