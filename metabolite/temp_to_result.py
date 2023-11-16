import pickle
import typer

def main():
    temp_cache = 'cache/matcher_temp.pickle'
    try:
        with open(temp_cache, 'rb') as f:
            temp = pickle.load(f)
    except FileNotFoundError:
        temp = {}
    
    results = {mol: [ms] for mols in temp.values() for mol, ms in mols.items()}
    
    with open('cache/matches.pickle', 'wb+') as f:
        pickle.dump(results, f)

if __name__ == '__main__':
    typer.run(main)