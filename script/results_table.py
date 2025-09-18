# Convert the result files (pickle) into a one human-readable table.


import sys
import os
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../'))    # to the project root (linear_motif_cc2/)
sys.path.append(project_root)
from linear_motif_cc2.data.protein_database_1_class import Protein    

import pickle
import pandas as pd
import numpy as np


# Toggle filter to False if a full table of all predictions (including bad motifs) is needed
filter = True


# initialize a list that stores all data of the table, which will be later converted into a pd df
data_rows = []

# read results
results_path = 'results_20250713/'
results = os.listdir(results_path)

# read protein database
protein_database_path = '../../data/protein_database_1.pickle'
with open(protein_database_path, 'rb') as f:
    protein_database = pickle.load(f)

n=0
for result in results:
    n+=1
    print(f'\r{n}', end='')
    path = results_path + result

    with open(path, 'rb') as f:
        data = pickle.load(f)
    for run in data:
        # filer away unconverged runs
        if run['completed_iterations'] == 1000 and filter:
            continue

        for i, (protein_id, motif) in enumerate(run['motifs'].items()):
            # filter away dropped out motifs
            if run['weights'][i] == 0 and filter:
                continue
            
            start = run['positions'][i]
            length = len(motif)
            end = start + length
            info_content = run['avg_score']

            # features
            protein = protein_database[protein_id]
            plm = np.mean(protein.plm[start:end])
            disorder = np.mean(protein.disorder[start:end])
            solvent_acc = np.mean(protein.rsa[start:end])
            conservation = np.mean(protein.rlc[start:end])
            
            start += 1  # to convert from 0-based python indices to 1-based both inclusive indices, increase start by 1, keep end the same
            data_rows.append([protein_id, motif, start, end, length, info_content, plm, disorder, solvent_acc, conservation])

df = pd.DataFrame(data_rows, columns=['protein', 'motif', 'start', 'end', 'length', 'info_content', 'plm', 'disorder', 'solvent_acc', 'conservation'])

# sort based on proteins
df = df.sort_values('protein')



comments = '''# Result Table of Linear Motif Prediction
# This table listed the predicted motifs using our linear motif prediction tool, specifing the infomation content and the four feature scores, calculated by averaging all residues in the motif.  
# Note 1: The four feature scores are z-normalize.
# Note 2: The information content is calculated from the match states of the HMM, the calculation is described in the paper
# Note 3: Only the confident predictions are presented here. To see all the predictions, including the filtered ones, either toggle the filter variable to False in the python script and generate a new table, or directly access the original pickle files.
'''

with open('results_table_20250713.csv', 'w') as f:
    f.write(comments)
    df.to_csv(f, index=False)


df = pd.read_csv('results_table_20250713.csv', comment="#")
print(df)




# with open('results_20250713/A0AV96.pickle', 'rb') as f:
#     data = pickle.load(f)

# data = data[0]
# for k, v in data.items():
#     print(f'{k}\t{v}')
