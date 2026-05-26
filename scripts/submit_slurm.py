import subprocess
import pickle
import pandas as pd
from io import StringIO
import os
import argparse


def submit(df, center, results_path):

    # estimate time
    size = df.at[center, 'size']
    if size < 50:
        estimated_time = '00:30:00'
    elif size < 100:
        estimated_time = '01:00:00'
    else:
        estimated_time = '02:00:00'

    # create bash file    
    bash_script = f"""#!/bin/bash
#SBATCH --job-name={center}          # Job name
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time={estimated_time}

module load python

echo "{center}"
python search.py {center} {results_path}/{center}.pickle
"""
    # Save to file
    bash_name = f'{center}.sh'
    with open(bash_name, 'w') as f:
        f.write(bash_script)

    # submit
    result = subprocess.run(['sbatch', bash_name], capture_output=True, text=True)
    os.remove(bash_name)    # delete
    output = result.stdout
    job_id = int(output.split()[-1])
    df.at[center, 'job_id'] = job_id

    return df

def initialize(biogrid_path):
    with open(biogrid_path, 'rb') as f:
        biogrid = pickle.load(f)
    
    # initialize
    df = pd.DataFrame(columns=['job_id', 'submitted', 'in_sq', 'size', 'result'])
    # populate
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

        df.loc[center] = {
            'job_id': pd.NA,
            'submitted': False,
            'in_sq': False,
            'size': size,
            'result': False
        }

    # sort
    df = df.sort_values('size', ascending=False)
    print(f'Initialized. Number of networks: {len(df)}')

    return df


def check(df, results_path):
    # check results: set all results as False, then make True for results
    df['result'] = False 
    
    results = os.listdir(results_path)
    results = [result.split('.')[0] for result in results]
    print(f'Completed results files: {len(results)}')

    for result in results:
        df.at[result, 'result'] = True

    # check queue
    df['in_sq'] = False
    queue = subprocess.run(['squeue'], capture_output=True, text=True)
    df_q = pd.read_fwf(StringIO(queue.stdout))
    in_q = list(df_q['NAME'])

    n=0
    for center, row in df.iterrows():
        if center in in_q: 
            df.at[center, 'in_sq'] = True        
            n+=1
    print(f'Jobs in sq currently: {n}')

    return df


def refill(df, max_jobs, results_path):
    # count how many in sq now
    n_in_sq = sum(df['in_sq'])
    print(f'Jobs in sq currently: {n_in_sq}')
    print(f'Target max job submission: {max_jobs}')

    df = df.sort_values('size', ascending=False)    # sort again just in case
    n = 0
    for center, row in df.iterrows():
        if n_in_sq >= max_jobs:
            break
        if not row['submitted'] and not row['result']:
            try:
                df = submit(df=df, center=center, results_path=results_path)
                df.at[center, 'submitted'] = True
                df.at[center, 'in_sq'] = True        
                n += 1
                n_in_sq += 1
            except Exception as e:
                print(f'Submission error at {center}')
                print(e)
                print()
    print(f'Jobs submitted: {n}')
    
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Auto-submission for SLURM High Performance Computing system",
        usage="python submit_slurm.py <biogrid_path> <results_folder> <dataframe_path>",
        epilog="Example: python submit_slurm.py ../interactome/biogrid_net.gpickle results_yyyymmdd/ df.pickle"
    )
    parser.add_argument("biogrid_path", type=str,
                        help="Path to the BioGRID network file in .gpickle format")
    parser.add_argument("results_folder", type=str,
                        help="Directory path where results will be saved")
    parser.add_argument("df_path", type=str,
                        help="Path to the pandas dataframe pickle file, which is used to keep track of the running status")
    args = parser.parse_args()

    df_path = args.df_path
    biogrid_path = args.biogrid_path
    
    if os.path.exists(df_path):
        with open(df_path, 'rb') as f:
            df = pickle.load(f)
        print("Loaded existing dataframe")
    else:
        df = initialize(biogrid_path=biogrid_path)
        with open(df_path, 'wb') as f:
            pickle.dump(df, f)
        print("Initialized new dataframe")

    results_path = args.results_folder
    os.makedirs(results_path, exist_ok=True)
    df = check(df=df, results_path=results_path)
    df = refill(df=df, max_jobs=950, results_path=results_path)
    print(df)

    with open(df_path, 'wb') as f:
        pickle.dump(df, f)


if __name__ == '__main__':
    main()
