
from functools import partial
from itertools import combinations
import yaml
from pathlib import Path
from typing import Dict, List
import pickle
from itertools import islice

import typer
import datamol as dm
import typer
from mpire import WorkerPool
from rdkit.Chem import Descriptors, DataStructs, rdFMCS

from metabolite.model import Metabolite

def filter_matches_ring_frac(matches, max_frac=0.8):
    results = []
    for m, mol, fp in matches:
        ring_info = mol.GetRingInfo()
        ring_atoms = set()
        for r in ring_info.AtomRings():
            ring_atoms.update(r)
        if len(ring_atoms) <= mol.GetNumAtoms() * 0.8:
            results.append((m, mol, fp))
    return results

def filter_matches_ro3(matches):
    results = []
    for m, mol, fp in matches:
        mw = Descriptors.ExactMolWt(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rb = Descriptors.NumRotatableBonds(mol)
        logp = Descriptors.MolLogP(mol)
        psa = Descriptors.TPSA(mol)
        if mw <= 300 and hbd <= 3 and hba <= 3 and rb <= 3 and logp <= 3 and psa <= 60:
            results.append((m, mol, fp))
    return results

def filter_matches_ro3(matches):
    results = []
    total_hbd = 0
    total_hba = 0
    for m, mol, fp in matches:
        mw = Descriptors.ExactMolWt(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rb = Descriptors.NumRotatableBonds(mol)
        logp = Descriptors.MolLogP(mol)
        psa = Descriptors.TPSA(mol)
        if mw <= 300 and hbd <= 3 and hba <= 3 and rb <= 3 and logp <= 3 and psa <= 60:
            results.append((m, mol, fp))
            total_hbd += hbd
            total_hba += hba
    return results

def squash_matches(matches, max_similarity=0.7, max_substructure=0.9):
    pairs = matches.copy()
    
    while len(pairs) > 1:
        for (a, a_mol, a_fp), (b, b_mol, b_fp) in combinations(pairs, 2):
            a_atoms = a_mol.GetNumAtoms()
            b_atoms = b_mol.GetNumAtoms()
            tan = DataStructs.TanimotoSimilarity(a_fp, b_fp)
            if tan > max_similarity:
                if a_atoms > b_atoms:
                    pairs.remove((a, a_mol, a_fp)) 
                else:
                    pairs.remove((b, b_mol, b_fp)) 
                break
            mcs = rdFMCS.FindMCS([a_mol, b_mol], ringMatchesRingOnly=True, completeRingsOnly=True, matchChiralTag=True)
            if mcs.numAtoms > min(a_atoms, b_atoms) * max_substructure:
                #print(a.smiles, b.smiles, mcs.numAtoms, min(a_atoms, b_atoms))
                if a_atoms > b_atoms:
                    pairs.remove((b, b_mol, b_fp))
                else:
                    pairs.remove((a, a_mol, a_fp))
                break
        else:
            break
    return pairs

def frag_fraction(mol, frags):
    matched = set()
    for f in frags:
        matched.update(mol.GetSubstructMatch(f))
    return len(matched) / mol.GetNumAtoms()

def is_druglike(mol, max_mw=700, max_rb=11):
    mol_mw = Descriptors.ExactMolWt(mol)
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    if mol_mw > max_mw or rotatable_bonds > max_rb:
        return False
    return True

def chunks(data, size=10000):
   it = iter(data)
   for i in range(0, len(data), size):
      yield {k: data[k] for k in islice(it, size)}

def refine_func(mol_matches: Dict[str, List[List[Metabolite]]], max_mw=700, max_rb=11, min_frac=0.7, max_similarity=0.7, max_substructure=0.9, max_ring_frac=0.8, ro3=True):
    f_mol_matches = {}
    for smiles, matches_list in mol_matches.items():
        mol = dm.to_mol(smiles)
        if not is_druglike(mol, max_mw=max_mw, max_rb=max_rb):
            #print('Filtered: druglike')
            continue
        f_matches_list = []
        for matches in matches_list:
            matches = filter_matches_ring_frac(matches, max_frac=max_ring_frac)
            if ro3:
                matches = filter_matches_ro3(matches)
            matches = squash_matches(matches, max_similarity=max_similarity, max_substructure=max_substructure)
            if len(matches) < 2:
                #print('Filtered: low matches after squash')
                continue
            m_mols = [m[1] for m in matches]
            ff = frag_fraction(mol, m_mols)
            if ff < min_frac:
                #print('Filtered: low frac')
                continue
            f_matches_list.append((matches, ff))
        if f_matches_list:
            f_matches = next(iter(sorted(f_matches_list, reverse=True, key=lambda x: x[1])))[0]
            f_mol_matches[smiles] = f_matches
    return f_mol_matches

def main(max_mw:int=600, max_rb:int=11, min_frac:float=0.7, max_similarity:float=0.7, max_substructure:float=0.9, max_ring_frac:float=0.8, ro3=True, cpus: int = None): 
    with open(f'cache/matches.pickle', 'rb') as f:
        mol_matches = pickle.load(f)
    
    Path('cache').mkdir(exist_ok=True)
    spec = f'refine_{max_mw}_{max_rb}_{min_frac}_{max_similarity}_{max_substructure}_{max_ring_frac}'
    if ro3:
        spec += f'_ro3'
    temp_cache = f'cache/{spec}.pickle'
    try:
        with open(temp_cache, 'rb') as f:
            temp = pickle.load(f)
    except FileNotFoundError:
        temp = {}
        
    print('Loaded')

    func = partial(refine_func, max_mw=max_mw, max_rb=max_rb, min_frac=min_frac, max_similarity=max_similarity, max_substructure=max_substructure, max_ring_frac=max_ring_frac, ro3=ro3)

    with WorkerPool(n_jobs=cpus) as pool:
        for mms in pool.imap_unordered(func, [(d,) for d in chunks(mol_matches)], progress_bar=True, progress_bar_options={'desc': 'Refining...'}):
            temp.update(mms)
            with open(temp_cache, 'wb+') as f:
                pickle.dump(temp, f)
    
    if len(temp):
        with open(temp_cache, 'wb+') as f:
            pickle.dump(temp, f)
        
        output = {}
        for smiles, matches in temp.items():
            output[smiles] = [[str(m.hmdb_id), str(m.smiles), 'Pathways: ' + '; '.join([f'{p.name} ({p.kegg_id})' if p.kegg_id else str(p.name) for p in m.pathways])] for m, _, _ in matches]
        
        with open(f'{spec}.yaml', 'w+') as f:
            yaml.dump(output, f)
        
    print(f'Result: {len(temp)} matches.')
    
if __name__ == '__main__':
    typer.run(main)