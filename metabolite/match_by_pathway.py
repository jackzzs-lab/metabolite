from itertools import combinations
from pathlib import Path
import pickle
import functools

from mpire import WorkerPool, cpu_count
import datamol as dm
from tqdm import tqdm
import typer
from rdkit import Chem
from rdkit.Chem import Descriptors, rdFMCS, DataStructs

dm.disable_rdkit_log()

fused_pattern = Chem.MolFromSmarts("[R&!R1]")
def fused_atoms(mol):
    return mol.GetSubstructMatches(fused_pattern)

def has_ring(mol):
    return len(Chem.GetSSSR(mol))

def is_valid_mcs(mol, m_mol):
    mcs = rdFMCS.FindMCS([mol, m_mol], ringMatchesRingOnly=True, completeRingsOnly=True, matchChiralTag=True)
    return mcs.numAtoms > int(m_mol.GetNumAtoms() * 0.9)

def get_ringsys(mol):
    ring_info = mol.GetRingInfo()
    ringsys = [set(ring) for ring in ring_info.AtomRings()]
    
    while True:
        for a, b in combinations(ringsys, 2):
            if len(a & b) >= 2:
                a.update(b)
                ringsys.remove(b)
                break
        else:
            break
    return ringsys

def squash_matches(matches):
    pairs = [(m, m_mol, Chem.RDKFingerprint(m_mol)) for m, m_mol in matches]
    
    while len(pairs) > 1:
        for (a, a_mol, a_fp), (b, b_mol, b_fp) in combinations(pairs, 2):
            tan = DataStructs.TanimotoSimilarity(a_fp, b_fp)
            if tan > 0.7:
                if a_mol.GetNumAtoms() < b_mol.GetNumAtoms():
                    pairs.remove((a, a_mol, a_fp)) 
                else:
                    pairs.remove((b, b_mol, b_fp)) 
                break
        else:
            break
    return pairs

def frags_matching_sub(db, metabolites_list):
    results = {}
    Path('runlog').mkdir(exist_ok=True)
    with open(f'runlog/{Path(db).stem}.run.log', 'w+', buffering=1) as runlog:
        with open(db, 'r') as f:
            n_matched = 0
            n_filtered = 0
            n_merged = 0
            for i, l in enumerate(f):
                try:
                    if not l.strip():
                        continue
                    mol_smiles = l.split()[0]
                    mol = dm.to_mol(mol_smiles)
                    if not mol:
                        continue
                    if not fused_atoms(mol):
                        continue
                    ringsys = get_ringsys(mol)
                    if any(len(rs) > 14 for rs in ringsys):
                        break
                    for metabolites in metabolites_list:
                        m_matches = []
                        for m, m_mol in metabolites:
                            try:
                                sub = set(mol.GetSubstructMatch(m_mol))
                                if sub:
                                    n_matched += 1
                                    for rs in ringsys:
                                        if sub & rs:
                                            if sub == rs:
                                                continue
                                            elif not rs.issubset(sub):
                                                n_filtered += 1
                                                break
                                    else:
                                        n_matched += 1
                                        m_matches.append((m, m_mol))
                            except:
                                continue
                        if len(m_matches) > 1:
                            for m, m_mol in m_matches:
                                if fused_atoms(m_mol):
                                    break
                            else:
                                continue
                            m_matches = squash_matches(m_matches)
                            if len(m_matches) > 3:
                                continue
                            elif len(m_matches) > 1:
                                if mol_smiles in results:
                                    results[mol_smiles].append(m_matches)
                                else:
                                    results[mol_smiles] = [m_matches]
                                m_smiles = [m[0].smiles for m in m_matches]
                                runlog.write(f'{mol_smiles} = {" + ".join(m_smiles)}\n')
                                n_merged += 1
                except Exception as e:
                    tqdm.write(f'Error: {e}')
                if i and i % 100 == 0:
                    tqdm.write(f'{db}: {i} processed, {n_matched} matched pair, {n_filtered} filtered pair, {n_merged} final merged.')
    return db, results

def main(database: Path, min_ms: int = 2, cpus: int = cpu_count()):
    database = list(Path(database).glob('*.smi'))
    print(f'Total: {len(database)} files.')
    
    Path('cache').mkdir(exist_ok=True)
    temp_cache = Path(f'cache/matcher_temp.pickle')
    try:
        with open(temp_cache, 'rb') as f:
            temp = pickle.load(f)
    except FileNotFoundError:
        temp = {}

    metabolites_list_cache = f'cache/matcher_metabolites_list_{min_ms}.pickle'
    try:
        with open(metabolites_list_cache, 'rb') as f:
            metabolites_list = pickle.load(f)
    except FileNotFoundError:
        with open('cache/hmdb_metabolites_pathways.pickle', 'rb') as f:
            pathways = pickle.load(f)
        if min_ms:
            pathways = {p: ms for p, ms in pathways.items() if len(ms) > min_ms}
        metabolites_list = []
        count = 0
        for ms in tqdm(list(pathways.values()), desc='Loaded pathways'):
            metabolites = []
            for m in ms:
                try:
                    m_mol = dm.to_mol(m.smiles)
                    if m_mol.GetNumAtoms() < 10:
                        continue
                    if not has_ring(m_mol):
                        continue
                    rotatable_bonds = Descriptors.NumRotatableBonds(m_mol)
                    if rotatable_bonds > 3:
                        continue
                except:
                    continue
                else:
                    metabolites.append((m, m_mol))
            metabolites_list.append(metabolites)
            count += len(metabolites)
        print(f'Total: {count} metabolites.')
        with open(metabolites_list_cache, 'wb+') as f:
            pickle.dump(metabolites_list, f)
    
    with WorkerPool(n_jobs=cpus) as pool:
        databases_todo = [db for db in database if db not in temp]
        func = functools.partial(frags_matching_sub, metabolites_list=metabolites_list)
        for db, mols in pool.imap_unordered(func, databases_todo, progress_bar=True, progress_bar_options={'desc': 'Matched database splits'}):
            temp[db] = mols
            tqdm.write(f'Finished: {db}')
            with open(temp_cache, 'wb+') as f:
                pickle.dump(temp, f)
    
    with open(temp_cache, 'wb+') as f:
        pickle.dump(temp, f)
    
    results = {mol: mss for mols in temp.values() for mol, mss in mols.items()}
    
    with open('cache/matches.pickle', 'wb+') as f:
        pickle.dump(results, f)
    
    print(f'Result: {len(results)} matches.')
    
if __name__ == '__main__':
    typer.run(main)
