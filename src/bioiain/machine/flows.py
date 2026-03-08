import os, sys, json, importlib, datetime
from ..utilities.logging import *
from .datasets import *
from .embeddings import *
from ..symmetries.elements import Monomer
from ..symmetries.interactions import PredictedMonomerContacts
from ..visualisation.pymol import PymolScript
from ..biopython.imports import loadPDB
import numpy as np
from ..machine import models



def predict(file_path, model_data_path, chain_id="A", use_temp=True, force=False, timestamp=False, foldseek_cmd="foldseek"):
    prediction_folder = "./predictions"
    if use_temp:
        prediction_folder = "/tmp/bioiain/predictions"
    seed = 6
    import random, torch

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # For multi-GPU


    name = os.path.basename(file_path).split(".")[0]
    if type(timestamp) is str:
        name = name + "_" + timestamp
    elif timestamp:
        name += datetime.datetime.now().strftime('_%y-%m-%d_%H-%M-%S')
    prediction_folder = os.path.join(prediction_folder, name)



    os.makedirs(prediction_folder, exist_ok=True)

    structure = loadPDB(file_path, name=name)
    print(structure)
    chain = None
    for ch in structure.get_chains():
        if ch.id == chain_id:
            chain = ch
            break
    assert chain is not None
    print(chain)
    monomer = Monomer.cast(chain)
    monomer.export(folder=prediction_folder)
    print(monomer)
    embedding = SaProtEmbedding(entity=monomer, folder=prediction_folder, force=force, foldseek_cmd=foldseek_cmd)
    data = json.load(open(model_data_path))
    print("Model class (data):", data.get("model"))
    model_class = getattr(models, data.get("model"))
    weights = data.get("weights", [])
    model = model_class(name="interactions", in_shape=data["in_shape"], num_classes=data["num_classes"],
                        weights=weights, inference=True)
    model.load(model_data_path, weights_only=True)
    # print(model.set_mode("no-dropout"))
    print(repr(model))
    # print(json.dumps(model.data, indent=4))

    pred_dataset = EmbeddingDataset(name=name, folder=prediction_folder)
    pred_dataset.add(embedding=embedding, key=monomer.get_name())
    label_to_index = model.data["label_to_index"]
    index_to_label = model.data["index_to_label"]
    print(pred_dataset)

    full_pred = []
    with torch.no_grad():
        for item in pred_dataset:
            out = model(item.t)
            print(out)
            # print(out.shape)
            if out.shape[0] > 2:
                pred = out.max(dim=0)[1]
            else:
                pred = out
            # print(pred)
            if len(label_to_index) > 2:
                p = index_to_label[str(pred.item())]
                full_pred.append(p)
            elif len(label_to_index) == 2:
                cp = round(pred[0].item() * 100)
                op = int(pred[1].item() > 0.5)
                full_pred.append(pred[1].item())
            else:
                full_pred.append(pred.item())
    print(full_pred)

    interaction = PredictedMonomerContacts(monomer, full_pred, label_to_index)
    # interaction._agglomerate()
    # interaction.plot_mpl()
    pred_path = interaction.save_structure(prediction_folder)
    script = PymolScript(name=f"{monomer.get_name()}_{chain.id}_prediction_pml_session",
                         folder=prediction_folder)
    script.load(pred_path, monomer.get_name())
    script.spectrum(monomer.get_name())
    if label_to_index is not None:
        script.print(json.dumps(label_to_index, indent=4))
    session_path = script.write_script()
    print("Session saved at:")
    print("pymol", session_path)
    # script.execute()
    return {"script": script, "prediction": pred_path, "timestamp":timestamp, "name":name,
            "label_to_index":label_to_index, "chain":chain, "session":session_path}









def display_monomer_labels(monomer_id, dataset, label_name, script=None, use_temp=True):


    folder = "./visualisations"
    if use_temp:
        folder = "/tmp/bioiain/visualisations"

    dataset.use_label(label_name)
    dataset.map()


    log("start", "VISUALISATION")
    log("header", f"Displaying monomer: {monomer_id}")

    embedding = dataset.embeddings[monomer_id]
    print(json.dumps(embedding, indent=4))

    label_path = embedding[label_name]
    with open(label_path, "r") as f:
        label = f.read().split(",")

    monomer = Monomer.recover(
        data_path=os.path.join(dataset.data["export_folder"], monomer_id.split("_")[0], "monomers", monomer_id))

    log(1, f"label: {label}")
    print(dataset.data["label_to_index"])
    if len(label) == 1:
        label = label[0]
    interaction = PredictedMonomerContacts(monomer, label, label_to_index=dataset.data["label_to_index"])
    view_path = interaction.save_structure(folder, extra_name="_visualised_monomer_contacts")
    log(1, interaction)

    if script is None:
        script = PymolScript(name=f"{monomer.get_name()}_visualisation_pml_session",
                             folder=folder)
    script.load(view_path, monomer.get_name())
    script.spectrum(monomer.get_name())

    session_path = script.write_script()
    print("Session saved at:")
    print("pymol", session_path)
    return script



