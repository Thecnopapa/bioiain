import os, json, sys

sys.path.append('..')


from src.bioiain.utilities import *
log("start", "aleph.py")

from src.bioiain.utilities.logging import *
tracemalloc_start()


from src.bioiain.aleph import *
from src.bioiain.base import *
from src.bioiain.machine import *
from src.bioiain.utilities.parallel import *

import torch, random
import  numpy as np

torch.set_num_threads(avail_cpus)
log(1, f"Torch using {avail_cpus} threads")

seed = 6
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)

if "monomers" in sys.argv:
    if not "--no-download" in sys.argv:
        DATA_FOLDER = downloadPDB("./data", "cath-monomeric",
                                  file_path="./data/cath-dataset-nonredundant-S20.monomeric.list",
                                  file_format="cif",
                                  overwrite=False)
    else:
        DATA_FOLDER = "./data/cath-monomeric"
    DATA_NAME = "monomers"
elif "receptors" in sys.argv:
    if not "--no-download" in sys.argv:
        DATA_FOLDER = downloadPDB("./data", "receptors",
                                  file_path="./data/receptors.txt",
                                  file_format="cif",
                                  overwrite=False)
    else:
        DATA_FOLDER = "./data/receptors"
    DATA_NAME = "receptors"
elif "lbds" in sys.argv:
    if not "--no-download" in sys.argv:
        DATA_FOLDER = downloadPDB("./data", "lbds",
                                  file_path="./data/LBDs.txt",
                                  file_format="cif",
                                  overwrite=False)
    else:
        DATA_FOLDER = "./data/lbds"
    DATA_NAME = "lbds"

elif "consensus" in sys.argv:
    if not "--no-download" in sys.argv:
        DATA_FOLDER = downloadPDB("./data", "consensus",
                                  file_path="./data/consensus.txt",
                                  file_format="cif",
                                  overwrite=False)
    else:
        DATA_FOLDER = "./data/consensus"
    DATA_NAME = "consensus"
else:
    DATA_FOLDER = downloadPDB( list_name="aleph", pdb_list=["1M2Z", "3HBB", "6F63", "5LXN", "3brf", "6e52", "7t2y", "3kg2", "2GEJ", "2bis"], data_dir="./data" )
    DATA_NAME = "aleph"


V0 = False
V1 = False
V2 = False
V3 = False
V4 = False

VA = False
VB = False
VC = False



if "--vB" in sys.argv:
    V0 = True
    VB = True
    DATA_NAME += "_vB"
    EMBEDDING_CLASS = CVEmbeddingVB

elif "--vC" in sys.argv:
    V0 = True
    VC = True
    DATA_NAME += "_vC"
    EMBEDDING_CLASS = CVEmbeddingVC

elif "--v1" in sys.argv:
    V1 = True
    VA = True
    DATA_NAME += "_v1"
    EMBEDDING_CLASS = CVEmbeddingV1
elif "--v2" in sys.argv:
    V2 = True
    VA = True
    DATA_NAME += "_v2"
    EMBEDDING_CLASS = CVEmbeddingV2

elif "--v1C" in sys.argv:
    V1 = True
    VC = True
    DATA_NAME += "_v1C"
    EMBEDDING_CLASS = CVEmbeddingV1C

elif"--v2C" in sys.argv:
    V2 = True
    VC = True
    DATA_NAME += "_v2C"
    EMBEDDING_CLASS = CVEmbeddingV2C

elif "--v3" in sys.argv or "--v3C" in sys.argv:
    V3 = True
    VC = True
    DATA_NAME += "_v3C"
    EMBEDDING_CLASS = CVEmbeddingV3C

elif "--v4" in sys.argv or "--v4C" in sys.argv:
    V4 = True
    VC = True
    DATA_NAME += "_v4C"
    EMBEDDING_CLASS = CVEmbeddingV4C

else:
    DATA_NAME += "_v0"
    EMBEDDING_CLASS = CVEmbedding
    V0 = True
    VA = True


log(1, "DATA NAME:", DATA_NAME)

if "-p" not in sys.argv:

    log("start", "Embeddings")
    log("title", "Embeddings")

    LR = 0.0001
    if "--lr" in sys.argv:
        LR = float(sys.argv[sys.argv.index("--lr") + 1])
    log(1, f"Learning rate: {LR}")

    MODEL_NAME = "Hope"
    if "--model" in sys.argv:
        MODEL_NAME = sys.argv[sys.argv.index("--model") + 1]
    MODEL_CLASS = getattr(models, MODEL_NAME)
    log(1, f"Model: {MODEL_CLASS}")


    dataset = EmbeddingDataset(name=f"tokens_{DATA_NAME}")

    if not ("--rebuild" in sys.argv or "--force" in sys.argv):
        dataset.load()
    log(2, dataset)
    total_files = len(os.listdir(DATA_FOLDER))
    if len(dataset) == 0:
        pool = None
        if "--thread" not in sys.argv:
            parts = [os.listdir(DATA_FOLDER)]
        else:
            parts = split_iterable(os.listdir(DATA_FOLDER), n_parts="half")
            pool = ThreadPool()

        def generate_embeddings(file_list=None):
            log("header", f"Generating embeddings... ({len(file_list)})")
            for n, file in enumerate(file_list):
                log("header", file)
                log("title", f"{dataset.n_ids()+1:4d}/{total_files:4d} ({file.split('.')[0]})")

                path = os.path.join(DATA_FOLDER, file)
                entity = BIEntity.from_file(path)

                if entity is None:
                    continue
                if len(entity) > 2000:
                    log("Warning", "entity too large!")
                    continue

                embedding = EMBEDDING_CLASS(entity=entity).embedding(force="--force" in sys.argv)

                if embedding is None:
                    log("warning", "No embedding for file:", file)
                    continue
                if "1M2Z" in file:
                    [print(r, e) for r, e in zip(entity.residues(), embedding.tensor())]
                    # entity.fragment().show_cvectors()

                print(embedding)
                dataset.add(embedding, key=entity.name())
                print(dataset)

        if len(parts) == 1 or pool is None:
            generate_embeddings(parts[0])
        else:
            for part in parts:
                pool.add(generate_embeddings, file_list=part)
            pool.start(wait=True)


        dataset.save()
    dataset.sequence_db()
    dataset.cluster(reassign=True, verbosity=3)

    dataset.save()
    dataset.align(verbose=True, build_tree=True)
    dataset.save()
    log("end", "Embeddings")


model = None
if "-t" in sys.argv:
    log("start", "Training")
    log("title", "Training")


    epochs = 100
    print(dataset)

    model = MODEL_CLASS(name=DATA_NAME, in_shape=dataset.get(0).t.shape, batch_size=0, lr=LR, embedding_class = EMBEDDING_CLASS)
    model.add_text("data", model.json())
    model.add_text("hparams", json.dumps({
        "model_name": model.__class__.__name__,
        "dataset": str(dataset),
        "label": dataset.data["label_key"],
        "seed": seed,
        "optimiser": model.optimisers.get(model.mode, "default").__class__.__name__,
        "loss_fn": model.criterions.get(model.mode, "default").__class__.__name__,
        "lr": LR,
        "batch_size": 0,
        "target_epochs": epochs,
        "device": DEVICE,
        }, indent=4))

    model.set_mode("autoencoder")
    model.mount()

    for n in range(epochs):
        log("start", "EPOCH", n)
        log("title", "EPOCH", n)
        model.set_mode("autoencoder")

        if "--no-plot" in sys.argv:
            model.plot_latent_space(dataset=None)
        else:
            model.plot_latent_space(dataset=dataset, max_points=1000, mesh_points=20)

        model.plot_tokens()



        n_items = len(dataset)
        for i, item in enumerate(dataset):

            loss, encoder_loss, decoder_loss = model(item.t)

            #loss = model.train(item, i, n_items)

            print(f"{i}/{n_items} LOSS: {loss.item():7.3f} ({encoder_loss.item():7.3f}/{decoder_loss.item():7.3f}) av:{model.running_loss[model.mode]/model.running_loss["total"]:7.3f}", end="\r")


            if (i+1) % 100000 == 0:
                logging.tracemalloc_top()



        #model.write_loss()
        #model.draw_all_tokens()
        model.save(temp=True)
        if not "--local" in sys.argv and ( (n+1) % 10 == 0 or (n+1)==epochs):
            model.send_run(host="iainvisa.com", key=os.environ.get("IAINVISA_FILE_KEY", None), epoch=n)
        model.add_epoch()
    model.save()
    log("end", "Training")



if "--tokenise" in sys.argv or "-t" in sys.argv:
    log("start", "Tokenisation")
    log("title", "Tokenisation")

    if model is None:
        log("header", "Loading saved model...")
        model_data_path = sys.argv[sys.argv.index("--md") + 1]
        log(1, "Model path:", model_data_path)
        data = json.load(open(model_data_path))
        log(1, "Model class (data):", data.get("model"))
        model_class = getattr(models, data.get("model"))

        model = model_class(name="inference", in_shape=data.get("in_shape"), inference=True)
        model.load(model_data_path)
    else:
        log("header", "Using loaded model...")

    log(1, "Model:", model)
    log(1, "Dataset:", dataset)

    tok_fasta = model._tokenise(dataset)
    matrix_path = model._build_blossum()
    model._align_tokens(dataset, tok_fasta, matrix="path", matrix_path=matrix_path, force=True)
    #model._align_tokens(dataset, tok_fasta)



if "-p" in sys.argv:
    log("start", "Prediction")
    log("title", "Prediction")
    with torch.no_grad():
        from src.bioiain.visualisation.pymol import PymolScript
        from src.bioiain.visualisation.plots import mpl_colours
        from src.bioiain.base.mmcif import write_atoms



        filepath = sys.argv[sys.argv.index("--file") + 1]
        model_data_path = sys.argv[sys.argv.index("--md") + 1]

        data = json.load(open(model_data_path))
        print("Model class (data):", data.get("model"))
        model_class = getattr(models, data.get("model"))

        model = model_class(name="inference", in_shape=data.get("in_shape"), inference=True)
        model.load(model_data_path)

        name = os.path.basename(filepath).split(".")[0]

        prediction_name = name +"_"+ datetime.datetime.now().strftime('_%y-%m-%d_%H-%M-%S')
        prediction_folder = os.path.join(SUBDIR_NAME, f"predictions/{prediction_name}")
        os.makedirs(prediction_folder, exist_ok=True)

        entity = BIEntity.from_file(filepath, code=name, force=True, export_folder=prediction_folder)

        entity.export()

        script = PymolScript(name=f"{name}", folder=prediction_folder)
        script.load(entity.path(minimal=False), entity.name())


        print(entity)
        print(len(entity.residues()))
        if model.data["embedding_class"] is not None:
            if EMBEDDING_CLASS.__name__ == model.data["embedding_class"]:
                log("warning", f"Embedding class ({EMBEDDING_CLASS.__name__}) does not match the model embedding class ({model.data['embedding_class']})")

        embedding = EMBEDDING_CLASS(entity=entity).embedding(force=True)
        print(embedding)
        
        entity = embedding.entity
        print(entity)
        residues = entity.residues()
        cvectors = entity.cvectors()

        print(len(residues), len(cvectors))



        entity.show_cvectors(script=script, execute=False)
        script.spectrum(entity.name(), color="blue_yellow_red")

        paths_emb = []
        t = embedding.tensor()
        for i in range(t.shape[-1]):
            for res in residues:
                res.set_bfactor(0)
            for cv, emb in zip(cvectors, t):
                res = cv.res2
                res.set_bfactor(emb[i])
            paths_emb.append(entity.export(sufix=f"EMB{i}"))
            script.load(paths_emb[-1])

        script.spectrum("*EMB*", color="blue_yellow_red", minimum=0, maximum=1)
        script.write_script()

        preds = [model._predict(e) for e in embedding.tensor()]

        print(len(cvectors), len(preds))
        assert len(cvectors) == len(preds)




        for res in residues:
            res.set_bfactor(0)


        tok_list = {}
        for cv, pred in zip(cvectors, preds): 
            res = cv.res2
            #print(res, pred[0])
            t = pred[0]
            res.set_bfactor(t)
            rname = f"TOK{t:2d}_{cv.chain}_{cv.resseq}"
            if not t in tok_list:
                tok_list[t] = [rname]
            else:
                tok_list[t].append(rname)
            cv_chain = BIChain.from_atoms([*cv.res1.atoms, *cv.res2.atoms, *cv.res3.atoms], code=rname, chain_id="A", complex=True, share=False)
            closest_chain = BIChain.from_atoms([*cv.closest.res1.atoms, *cv.closest.res2.atoms, *cv.closest.res3.atoms], code=rname, chain_id="B", complex=True, share=False)
            atoms = [*cv_chain.all_atoms(), *closest_chain.all_atoms()]
            script.load(write_atoms(atoms, os.path.join(TEMP_FOLDER, "trash", rname)), rname)

        path_tok = entity.export(sufix="tokens")
        script.load(path_tok)
        script.spectrum("*tokens", color="_".join(mpl_colours)+"_"+"_".join(mpl_colours), minimum=0, maximum=19)

        #script.disable("(all)")


        for tok, tok_res in tok_list.items():
            script.group(f"TOK{tok:2d}_", f"Token_{tok}")
            for res in tok_res:
                script.align(f"({res} and c. A)", f"({tok_res[0]} and c. A)")
        script.hide("TOK*")
        script.show("TOK*", "lines")
        script.spectrum("TOK*", color="_".join(mpl_colours) + "_" + "_".join(mpl_colours), minimum=0, maximum=19)

        script.write_script()

        model.plot_latent_space(plot_preds=preds, show=True, fig_dir=prediction_folder)

        script.execute(compile=True)
        script.execute()
    log("end", "Prediction")






log("end", "DONE")









