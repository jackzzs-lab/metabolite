import pickle
with open('cache/hmdb_metabolites_pathways.pickle', 'rb') as f:
    pathways = pickle.load(f)
    
from metabolite.model import Metabolite
results = {}

for i, (p, ms) in enumerate(pathways.items()):
    results[p] = [Metabolite(name=m.name, hmdb_id=m.hmdb_id, smiles=m.smiles, pathways=tuple(m.pathways), enzymes=tuple(m.enzymes), cas_registry_number=m.cas_registry_number) for m in ms]
    print(i)

with open('cache/hmdb_metabolites_pathways_2.pickle', 'wb+') as f:
    pickle.dump(results, f)