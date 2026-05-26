# Run entire proteome using user's own computer

import subprocess
import pickle
import argparse

parser = argparse.ArgumentParser(
    description="Run the linear motif searching algorithm on the entire proteome.",
    usage="python run_proteome.py <biogrid_path> <output_folder>",
    epilog="Example: python run_proteome.py ../interactome/biogrid_net.gpickle results_yyyymmdd/"
)
parser.add_argument("biogrid_path", type=str)
parser.add_argument("output_folder", type=str)
args = parser.parse_args()

biogrid_path = args.biogrid_path
results_folder = args.output_folder

# read biogrid
with open(biogrid_path, 'rb') as f:
    biogrid = pickle.load(f)

centers = []
for center in biogrid.nodes:
    if center == '-':
        continue

    # network size
    filtered = set()
    if len(biogrid[center]) <= 50 and len(biogrid[center]) >=10:
        for u, v, data in biogrid.edges(center, data=True):
            if v == '-' or u == '-':
                continue
            filtered.add(v)
    else:
        for u, v, data in biogrid.edges(center, data=True):
            if v == '-' or u == '-':
                continue
            if len(data['BiogridIDs']) >= 2 or data['low_throughput']:
                filtered.add(v)
    size = len(filtered)

    if size < 10 or size > 200:
        continue
    
    centers.append(center)


# run and save
for center in centers:
    try:
        subprocess.run(['python', 'search.py', f'{center}', f'{results_folder}{center}.pickle'])
    except Exception as e:
        print(f'Error occured at {center}\n{e}\n')
