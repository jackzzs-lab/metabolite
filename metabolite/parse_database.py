import pickle
from typing import List
from pathlib import Path

import typer
from lxml import etree
from tqdm import tqdm

from metabolite.model import Pathway, Enzyme, Metabolite

Def = object()
def get_xpath_first(element, xpath, default=Def, **kw):
    res = element.xpath(xpath, **kw)
    if res:
        return res[0]
    elif default is not Def:
        return default
    else:
        raise ValueError(f'XPath "{xpath}" not found in element {element.tag}')

def iter_metabolites():
    doc = etree.iterparse('data/hmdb_metabolites.xml', events=('start', 'end'))
    _, root = next(doc)
    start_tag = None
    for event, element in doc:
        if event == 'start' and start_tag is None:
            start_tag = element.tag
        if event == 'end' and element.tag == start_tag:
            yield element
            start_tag = None
            root.clear()

def parse_metabolite(element):
    smiles = get_xpath_first(element, "*[local-name()='smiles']/text()", None)
    if not smiles:
        return None
    pathways = []
    for e in element.xpath("*[local-name()='biological_properties']/*[local-name()='pathways']/*[local-name()='pathway']"):
        pathways.append(
            Pathway(
                name = get_xpath_first(e, "*[local-name()='name']/text()"),
                smpdb_id = get_xpath_first(e, "*[local-name()='smpdb_id']/text()", None),
                kegg_id = get_xpath_first(e, "*[local-name()='kegg_map_id']/text()", None),
            )
        )
    enzymes = []
    for e in element.xpath("*[local-name()='protein_associations']/*[local-name()='protein']"):
        uniprot_id = get_xpath_first(e, "*[local-name()='uniprot_id']/text()", None)
        if not uniprot_id:
            continue
        enzymes.append(
            Enzyme(
                name = get_xpath_first(e, "*[local-name()='name']/text()"),
                uniprot_id = uniprot_id,
                gene_name = get_xpath_first(e, "*[local-name()='gene_name']/text()", None),
            )
        )
    return Metabolite(
        name = get_xpath_first(element, "*[local-name()='name']/text()"),
        hmdb_id = get_xpath_first(element, "*[local-name()='accession']/text()"),
        smiles = smiles,
        cas_registry_number = get_xpath_first(element, "*[local-name()='cas_registry_number']/text()", None),
        pathways = pathways,
        enzymes=enzymes
    )

def main(): 
    hmdb_metabolites_cache = 'cache/hmdb_metabolites.pickle'
    if not Path(hmdb_metabolites_cache).exists():
        metabolites = []
        for element in tqdm(iter_metabolites(), desc='Loading...'):
            metabolite = parse_metabolite(element)
            if metabolite:
                metabolites.append(metabolite)
            
        with open(hmdb_metabolites_cache, 'wb+') as f:
            pickle.dump(metabolites, f)

    try:
        metabolites
    except NameError:
        with open(hmdb_metabolites_cache, 'rb') as f:
            metabolites: List[Metabolite] = pickle.load(f)

    pathways = {}
    for m in tqdm(metabolites, desc='Analyzing...'):
        for p in m.pathways:
            if p in pathways:
                pathways[p].append(m)
            else:
                pathways[p] = [m]
        
    with open('cache/hmdb_metabolites_pathways.pickle', 'wb+') as f:
            pickle.dump(pathways, f)

    print('Top pathways:')
    for i, j in enumerate(sorted(pathways, key=lambda k: len(pathways[k]), reverse=True)):
        if i == 0:
            top = len(pathways[j])
        print(f'  {i+1} {j.name} ({len(pathways[j])})')
        if i > 8:
            break

    lg = lambda k: len(list(p for p in pathways if len(pathways[p]) > k))
    lgstr = lambda k:f'>{top*k:.0f} ({k*100:.0f}%): {lg(top*k)}'
    print('\nStatistics:\n  ' + ', '.join(lgstr(k / 10) for k in range(1, 10)))
    
if __name__ == '__main__':
    typer.run(main)