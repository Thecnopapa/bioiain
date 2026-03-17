import os, json, sys

sys.path.append('..')


from src.bioiain.utilities import *
log("start", "aleph.py")




from src.bioiain.aleph import *
from src.bioiain.base import *
from src.bioiain.machine import *



from src.bioiain.biopython.imports import *
folder = downloadPDB( list_name="aleph", pdb_list=["1M2Z", "3HBB", "6F63", "5LXN", "3brf"], data_dir="./data" )



dataset = EmbeddingDataset(name="aleph_test")

for file in os.listdir(folder):

    log("header", file)
    log("title", file)

    path = os.path.join(folder, file)
    entity = BIEntity.from_file(path)

    embedding = CVEmbedding(entity=entity).embedding(force="--force" in sys.argv)
    print(embedding)
    dataset.add(embedding, key=entity.name())
    print(dataset)


dataset.save()
epochs = 100
model = Despair(name="test", in_shape=dataset[0].t.shape, dry=True)


for n in range(len(dataset)):
    item = dataset[n]
    print(item)
    out = model(item.t)
    print(out)











