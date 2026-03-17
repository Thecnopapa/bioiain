import os, json, sys

sys.path.append('..')


from src.bioiain.utilities import *
log("start", "aleph.py")




from src.bioiain.aleph import *
from src.bioiain.base import *
from src.bioiain.machine import *
from src.bioiain.utilities.parallel import *



from src.bioiain.biopython.imports import *
if "monomers" in sys.argv:
    if not "--no-download" in sys.argv:
        DATA_FOLDER = downloadPDB("./data", "cath-monomeric",
                                  file_path="./data/cath-dataset-nonredundant-S20.monomeric.list",
                                  file_format="cif",
                                  overwrite=False)
    else:
        DATA_FOLDER = "./data/cath-monomeric"
    DATA_NAME = "monomers"
else:
    DATA_FOLDER = downloadPDB( list_name="aleph", pdb_list=["1M2Z", "3HBB", "6F63", "5LXN", "3brf", "6e52", "7t2y"], data_dir="./data" )
    DATA_NAME = "aleph"

LR = 0.0001
if "--lr" in sys.argv:
    LR = float(sys.argv[sys.argv.index("--lr") + 1])

log(1, f"Learning rate: {LR}")

dataset = EmbeddingDataset(name="aleph_test")

if not "--rebuild" in sys.argv and not "--force" in sys.argv:
    dataset.load()
else:
    parts = split_iterable(os.listdir(DATA_FOLDER))
    pool = ThreadPool()

    def generate_embeddings(file_list=None):
        log("header", f"Generating embeddings... ({len(file_list)})")
        for file in file_list:
            log("header", file)
            log("title", file)

            path = os.path.join(DATA_FOLDER, file)
            entity = BIEntity.from_file(path)

            embedding = CVEmbedding(entity=entity).embedding(force="--force" in sys.argv)

            if embedding is None:
                log("warning", "No embedding for file:", file)
                continue
            if "7T2Y" in file:
                [print(r, e) for r, e in zip(entity.residues(), embedding.tensor())]
                # entity.fragment().show_cvectors()

            print(embedding)
            dataset.add(embedding, key=entity.name())
            print(dataset)


    for part in parts:
        pool.add(generate_embeddings, file_list=part)
    pool.start(wait=True)


    dataset.save()


if "-t" in sys.argv:

    epochs = 100
    print(dataset)

    model = Despair(name=DATA_NAME, in_shape=dataset.get(0).t.shape, batch_size=1, lr=LR)
    model.add_text("data", model.json())

    for n in range(epochs):
        log("start", "EPOCH", n)
        log("title", "EPOCH", n)
        model.cluster_latent_space(dataset)
        model.plot_current_state(dataset=dataset)

        for item in dataset:
            latent = model(item.t, to_latent=True)
            new_latent = model.get_closest_latent(latent)
            out = model(new_latent, from_latent=True)
            loss = model.loss(out, item, zero_optims="autoencoder", step="autoencoder")
            print(f"LOSS: {loss.item():7.3f} {model.running_loss["default"]/model.running_loss["total"]:7.3f}", end="\r")

        model.write_loss()
        model.save()
        model.add_epoch()



if "-p" in sys.argv:
    with torch.no_grad():
        filepath = sys.argv[sys.argv.index("--file") + 1]
        modelpath = sys.argv[sys.argv.index("--model") + 1]

        model = Despair(name="test", in_shape=dataset.get(0).t.shape, inference=True)
        model.load(modelpath)
        entity = BIEntity.from_file(filepath, code="pred", force=True)
        embedding = CVEmbedding(entity=entity).embedding(force=True)
        out = model(embedding.tensor(), to_latent=True)
        preds = model.predict_tokens(out)
        print(preds)





log("end", "DONE")









