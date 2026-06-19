import os, json, sys


sys.path.append('..')


from src.bioiain.utilities import *
log("start", "aleph.py")

from src.bioiain.utilities.logging import *
tracemalloc_start()


from src.bioiain.base import *
from src.bioiain.aleph import *
from src.bioiain.machine import *
from src.bioiain.utilities.parallel import *
import vqvae_models as models
import vqvae_embeddings as embeddings
from relative_contactability import calculate_relative_contactability

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
        DATA_FOLDER = downloadPDBlist("./data", "cath-monomeric",
                                      file_path="./data/cath-dataset-nonredundant-S20.monomeric.list",
                                      file_format="cif",
                                      overwrite=False)
    else:
        DATA_FOLDER = "./data/cath-monomeric"
    DATA_NAME = "monomers"
elif "receptors" in sys.argv:
    if not "--no-download" in sys.argv:
        DATA_FOLDER = downloadPDBlist("./data", "receptors",
                                      file_path="./data/receptors.txt",
                                      file_format="cif",
                                      overwrite=False)
    else:
        DATA_FOLDER = "./data/receptors"
    DATA_NAME = "receptors"
elif "lbds" in sys.argv:
    if not "--no-download" in sys.argv:
        DATA_FOLDER = downloadPDBlist("./data", "lbds",
                                      file_path="./data/LBDs.txt",
                                      file_format="cif",
                                      overwrite=False)
    else:
        DATA_FOLDER = "./data/lbds"
    DATA_NAME = "lbds"

elif "consensus" in sys.argv:
    if not "--no-download" in sys.argv:
        DATA_FOLDER = downloadPDBlist("./data", "consensus",
                                      file_path="./data/consensus.txt",
                                      file_format="cif",
                                      overwrite=False)
    else:
        DATA_FOLDER = "./data/consensus"
    DATA_NAME = "consensus"
else:
    DATA_FOLDER = downloadPDBlist(list_name="aleph", pdb_list=["1M2Z", "3HBB", "6F63", "5LXN", "3brf", "6e52", "7t2y", "3kg2", "2GEJ", "2bis"], data_dir="./data")
    DATA_NAME = "aleph"


FORCE = "--force" in sys.argv or "-f" in sys.argv
REBUILD = "--rebuild" in sys.argv or "-r" in sys.argv
RELATIVE = "--relative" in sys.argv or "--rel" in sys.argv
if RELATIVE:
    EMBEDDING_CLASS = embeddings.ExpandedALEPHEmbedding0
else:
    EMBEDDING_CLASS = embeddings.ALEPHEmbedding
DATA_NAME += "_"+ EMBEDDING_CLASS.__name__


log(1, "DATA NAME:", DATA_NAME)

if "-p" not in sys.argv:

    log("start", "Embeddings")
    log("title", "Embeddings")

    LR = 0.0001
    if "--lr" in sys.argv:
        LR = float(sys.argv[sys.argv.index("--lr") + 1])
    log(1, f"Learning rate: {LR}")

    MODEL_NAME = "Summer"
    if "--model" in sys.argv:
        MODEL_NAME = sys.argv[sys.argv.index("--model") + 1]
    MODEL_CLASS = getattr(models, MODEL_NAME)
    log(1, f"Model: {MODEL_CLASS}")

    DATASET_NAME = f"DATASET_{DATA_NAME}"
    dataset = EmbeddingDataset(name=DATASET_NAME)

    if not (REBUILD or FORCE):
        dataset.load()
    log(2, dataset)


    BLACKLIST = []
    blacklist_file = os.path.join(DATA_FOLDER, "BLACK.list")
    if os.path.exists(blacklist_file) and not FORCE:
        with open(blacklist_file, "r") as bl:
            for line in bl:
                BLACKLIST.append(line.strip().replace("\n", "").split(":")[0].strip())
    else:
        with open(blacklist_file, "w") as bl:
            bl.write("BLACK.list\n")


    file_list = os.listdir(DATA_FOLDER)
    total_files = len(file_list) - 1
    if not FORCE:
        log(1, "BLACKLIST:", BLACKLIST)
        file_list = [fl for fl in file_list if fl not in BLACKLIST]
        total_files = len(file_list)

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
                if file == "BLACK.list":
                    continue

                log("header", f"{dataset.n_ids()+1:4d}/{total_files:4d} ({file.split('.')[0]}) ({EMBEDDING_CLASS.__name__})")
                log("title", f"{dataset.n_ids()+1:3d}/{total_files:3d} ({EMBEDDING_CLASS.__name__})")

                if file in BLACKLIST and not FORCE:
                    log("warning", f"File in blacklist: {file}")
                    continue


                path = os.path.join(DATA_FOLDER, file)
                try:
                    entity = FragmentedStructure.from_file(path, no_atoms=True)
                except Exception as e:
                    log("Warning", "Skipping embedding for:", file, f"({e})")
                    with open(blacklist_file, "a") as f:
                        f.write(f"{file}:Entity load error\n")
                    continue
                log(1, "Pre-loading embedding...")
                embedding = embeddings.ALEPHProteinEmbedding(entity=entity, residue_embedding_class=EMBEDDING_CLASS, dry=True)
                log(2, embedding, f"EXISTS={embedding.exists()}")

                if not embedding.exists() or FORCE:
                    log(1, "Generating embedding...")

                    try:
                        entity = FragmentedStructure.from_file(path)
                    except Exception as e:
                        log("Warning", "Skipping embedding for:", file, f"({e})")
                        with open(blacklist_file, "a") as f:
                            f.write(f"{file}:Entity load error\n")
                    if len(entity) > 2000:
                        log("Warning", "Entity too large!")
                        with open(blacklist_file, "a") as f:
                            f.write(f"{file}:Too large\n")
                        continue
                    try:
                        entity.db(force=True)
                    except SequenceNotFound:
                        log("Warning", "Entity has no sequence!")
                        with open(blacklist_file, "a") as f:
                            f.write(f"{file}:No sequence\n")
                        continue

                    embedding = embeddings.ALEPHProteinEmbedding(entity=entity, residue_embedding_class=EMBEDDING_CLASS)

                    try:
                        embedding.generate()
                        embedding.save()
                        entity.export()
                    except (ALEPHError, NoEmbeddingForThisProtein, StructureLoadException) as e:
                        print(e)
                        embedding = None


                else:
                    log(1, "Embedding already generated")
                    embedding = embedding.reload()
                log(2, embedding)
                #print("#####")



                if embedding is None:
                    log("warning", "No embedding for file:", file)
                    with open(blacklist_file, "a") as f:
                        f.write(f"{file}:No embedding\n")
                    continue
                if "1M2Z" in file:
                    [print(r, e) for r, e in zip(entity.residues(), embedding.tensor())]
                    # entity.fragment().show_cvectors()

                #print(embedding)
                dataset.add(embedding, key=entity.name())
                dataset.save(temp=True)
                log(2, dataset)
                if (n+1) % 100 == 0:
                    tracemalloc_top()

        if len(parts) == 1 or pool is None:
            generate_embeddings(parts[0])
        else:
            for part in parts:
                pool.add(generate_embeddings, file_list=part)
            pool.start(wait=True)


        dataset.save()
        dataset.sequence_db(force=True)
        dataset.cluster(reassign=True, force=True)
        #dataset.align(verbose=True, build_tree=True, force=True)
        dataset.save()


    dataset.sequence_db()
    dataset.cluster(reassign=True)
    dataset.save()
    log("end", "Embeddings")


    if RELATIVE:
        log("start", "Relative Embeddings")

        DATASET_NAME = DATASET_NAME+"_relative"
        EMBEDDING_CLASS = embeddings.RelativeALEPHEmbedding
        relative = EmbeddingDataset(name=DATASET_NAME)
        log(2, relative)

        if not (REBUILD or FORCE):
            relative.load()

        if len(relative) == 0:

            relative = calculate_relative_contactability(dataset)
            relative.save()
            relative.sequence_db(force=True)
            relative.cluster(reassign=True, force=True)
            relative.save()
            log(2, relative)
        log("end", "Relative Embeddings")
        dataset = relative


model = None
if "-t" in sys.argv and not ("-p" in sys.argv):
    log("start", "Training")
    log("title", "Training")


    epochs = 50
    if "--epochs" in sys.argv:
        epochs = int(sys.argv[sys.argv.index("--epochs") + 1])

    n_dots = 500
    if "--n-dots" in sys.argv:
        n_dots = int(sys.argv[sys.argv.index("--n-dots") + 1])
    n_squares = 15
    if "--n-squares" in sys.argv:
        n_squares = int(sys.argv[sys.argv.index("--n-squares") + 1])

    log(1, "EMBEDDING_CLASS:", EMBEDDING_CLASS)
    log(1, "DATASET:", dataset)

    try:
        assert RELATIVE
        for k, e in dataset.embeddings.items():
            print(k, set([float(dataset[i].t[8].item()) for i in range(e["start"], e["end"])]))
    except:
        pass


    model = MODEL_CLASS(name=DATASET_NAME, in_shape=dataset.get(0).t.shape, batch_size=0, lr=LR, embedding_class = EMBEDDING_CLASS)
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
        "n_dots": n_dots,
        "n_squares": n_squares,
        }, indent=4))

    model.set_mode("autoencoder")
    model.mount()

    total_params = sum(p.numel() for p in model.submodels["autoencoder"].parameters())
    log(1, "Number of parameters in the model:", total_params)

    for n in range(epochs):
        log("start", "EPOCH", n, model.__class__.__name__, DATASET_NAME)
        log("title", "EPOCH", n, model.__class__.__name__, DATASET_NAME)
        model.set_mode("autoencoder")

        if "--no-plot" in sys.argv:
            model.plot_latent_space(dataset=None)
        else:
            model.plot_latent_space(dataset=dataset, max_points=n_dots, mesh_points=n_squares)
            model.plot_latent_dimensions(dataset=dataset, max_points=n_dots, r_threshold=5)
            model.plot_latent_dimensions(dataset=dataset, max_points=n_dots, r_threshold=100 , plot_raw=False)

        model.plot_tokens()



        n_items = len(dataset)
        for i, item in enumerate(dataset):

            loss, encoder_loss, decoder_loss = model(item.t)

            #loss = model.train(item, i, n_items)

            print(f"{i}/{n_items} LOSS: {loss.item():7.3f} ({encoder_loss.item():7.3f}/{decoder_loss.item():7.3f}) av:{model.running_loss[model.mode]/model.running_loss['total']:7.3f}", end="\r")


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
        print(embeddings)
        if model.data["embedding_class"] is not None:
            if EMBEDDING_CLASS.__name__ != model.data["embedding_class"]:
                log("warning", f"Embedding class ({EMBEDDING_CLASS.__name__}) does not match the model embedding class ({model.data['embedding_class']})")

        embedding = EMBEDDING_CLASS(entity=entity).embedding(force=True)
        dd = EmbeddingDataset(name = prediction_name, folder=prediction_folder)
        dd.add(embedding)
        print(embedding)
        
        entity = embedding.entity
        print(entity)
        residues = entity.residues()
        cvectors = entity.cvectors()

        print(len(residues), len(cvectors))



        entity.show_cvectors(script=script, execute=False)
        script.spectrum(entity.name(), color="blue_yellow_red")

        paths_emb = []
        paths_dec = []
        paths_disc = []
        t = embedding.tensor()
        preds = [model._predict(e) for e in embedding.tensor()]
        [print(p) for p in preds]
        decoded = [model._decode(pred[3].detach()) for pred in preds]
        discretised = [model._decode(pred[2].detach()) for pred in preds]

        names = ["len i", "len j", "angle ij", "dist ij", "dist lig", "contactability", "SASA", "dihedral",
                 "t1", "t2"]


        for i in range(t.shape[-1]):
            for res in residues:
                res.set_bfactor(0)
            for cv, emb in zip(cvectors, t):
                res = cv.res2
                res.set_bfactor(emb[i])
            paths_emb.append(entity.export(sufix=f"EMB_{names[i]}"))
            script.load(paths_emb[-1])

            for res in residues:
                res.set_bfactor(0)
            for cv, dec in zip(cvectors, decoded):
                res = cv.res2
                #print(dec)
                res.set_bfactor(dec[i])
            paths_dec.append(entity.export(sufix=f"DEC_{names[i]}"))
            script.load(paths_dec[-1])

            for res in residues:
                res.set_bfactor(0)
            for cv, dis in zip(cvectors, discretised):
                res = cv.res2
                #print(dis)
                res.set_bfactor(dis[0][i])
            paths_disc.append(entity.export(sufix=f"DIS_{names[i]}"))
            script.load(paths_disc[-1])


        script.spectrum("*EMB*", color="blue_yellow_red", minimum=0, maximum=1)
        script.spectrum("*DEC*", color="blue_yellow_red", minimum=0, maximum=1)
        script.spectrum("*DIS*", color="blue_yellow_red", minimum=0, maximum=1)
        for name in names:
            script.group(name, name, also_suffix=True)

        script.write_script()



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
            # rname = f"TOK{t:2d}_{cv.chain}_{cv.resseq}"
            # if not t in tok_list:
            #     tok_list[t] = [rname]
            # else:
            #     tok_list[t].append(rname)
            # cv_chain = BIChain.from_atoms([*cv.res1.atoms, *cv.res2.atoms, *cv.res3.atoms], code=rname, chain_id="A", complex=True, share=False)
            # closest_chain = BIChain.from_atoms([*cv.closest.res1.atoms, *cv.closest.res2.atoms, *cv.closest.res3.atoms], code=rname, chain_id="B", complex=True, share=False)
            # atoms = [*cv_chain.all_atoms(), *closest_chain.all_atoms()]
            # script.load(write_atoms(atoms, os.path.join(TEMP_FOLDER, "trash", rname)), rname)

        path_tok = entity.export(sufix="tokens")
        script.load(path_tok)
        script.spectrum("*tokens", color="_".join(mpl_colours)+"_"+"_".join(mpl_colours), minimum=0, maximum=19)


        # for tok, tok_res in tok_list.items():
        #     script.group(f"TOK{tok:2d}_", f"Token_{tok}")
        #     for res in tok_res:
        #         script.align(f"({res} and c. A)", f"({tok_res[0]} and c. A)")
        # script.hide("TOK*")
        # script.show("TOK*", "lines")
        # script.spectrum("TOK*", color="_".join(mpl_colours) + "_" + "_".join(mpl_colours), minimum=0, maximum=19)

        script.set("grid_mode", 1)

        script.write_script()

        model.plot_latent_space(dataset=dd, plot_preds=preds, show=True, fig_dir=prediction_folder, mesh_points=10)

        script.execute(compile=True)
        script.execute()
    log("end", "Prediction")






log("end", "DONE")









