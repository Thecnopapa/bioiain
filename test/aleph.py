import os, json, sys

sys.path.append('..')


from src.bioiain.utilities import *
log("start", "aleph.py")




from src.bioiain.aleph import *
from src.bioiain.base import *
from src.bioiain.machine import *



from src.bioiain.biopython.imports import *
folder = downloadPDB( list_name="aleph", pdb_list=["1M2Z", "3HBB", "6F63", "5LXN", "3brf", "6e52", "7t2y"], data_dir="./data" )



dataset = EmbeddingDataset(name="aleph_test")

for file in os.listdir(folder):

    log("header", file)
    log("title", file)

    path = os.path.join(folder, file)
    entity = BIEntity.from_file(path)

    embedding = CVEmbedding(entity=entity).embedding(force="--force" in sys.argv)
    if "7T2Y" in file:
        [print(r, e) for r, e in zip(entity.residues(), embedding.tensor())]
        #entity.fragment().show_cvectors()

    print(embedding)
    dataset.add(embedding, key=entity.name())
    print(dataset)


dataset.save()
epochs = 100
model = Despair(name="test", in_shape=dataset.get(0).t.shape)



for n in range(epochs):
    log("start", "EPOCH", n)
    log("title", "EPOCH", n)
    model.cluster_latent_space(dataset)
    model.plot_current_state(dataset=dataset)

    for item in dataset:
        latent = model(item.t, to_latent=True)
        new_latent = model.get_closest_latent(latent)
        out = model(new_latent, from_latent=True)

        print(f"LOSS: {model.loss(out, item).item():7.2f} {model.running_loss["default"]/model.running_loss["total"]:7.3f}", end="\r")


    model.plot_current_state(dataset=dataset)
    model.save()
    model.add_epoch()


log("end", "DONE")










