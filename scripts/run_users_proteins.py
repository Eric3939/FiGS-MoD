import numpy as np
import pickle
import argparse
import multiprocessing as mp
from functools import partial
from pathlib import Path
import pandas as pd
from figsmod.paths import PROTEIN_DB, RUNS, TABLES
from figsmod.search import load_protein_database, compute


def results_to_csv(results, input_proteins_id, output_csv_path):
    """Convert pickle results list to a CSV matching results_table.py format."""
    data_rows = []
    for run in results:
        if run['completed_iterations'] == 1000:
            continue
        for i, (protein_id, motif) in enumerate(run['motifs'].items()):
            if run['weights'][i] == 0:
                continue
            start = run['positions'][i]
            length = len(motif)
            end = start + length
            info_content = run['avg_score']
            data_rows.append([
                protein_id, motif, start + 1, end,
                length, round(info_content, 4),
                run['k'], run['seed']
            ])

    df = pd.DataFrame(data_rows, columns=[
        'protein', 'motif', 'start', 'end',
        'length', 'info_content', 'k', 'seed'
    ])
    df = df.sort_values('info_content', ascending=False)
    output_csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv_path, index=False)
    print(f"CSV results saved to {output_csv_path}")
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Run FiGS-MoD on a user-defined set of proteins.",
        usage="python run_users_proteins.py <proteins_file> [options]",
        epilog="Example: python run_users_proteins.py my_proteins.txt\n"
               "         python run_users_proteins.py my_proteins.txt --k 5 7 --runs 4",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "proteins_file", type=str,
        help="Path to a text file with one UniProt ID per line.\n"
             "These are the proteins hypothesised to share a common motif.\n"
             "Example contents:\n  Q15084\n  P12345\n  O95600"
    )
    parser.add_argument(
        "--output_dir", type=str, default=None,
        help="Directory to save outputs (default: runs/<proteins_file_stem>/)"
    )
    parser.add_argument(
        "--k", type=int, nargs="+", default=[7],
        help="Motif length(s) to search. Pass one or more values. Default: 7\n"
             "Example: --k 5 7 9"
    )
    parser.add_argument(
        "--runs", type=int, default=1,
        help="Number of independent runs per k value. Default: 1\n"
             "More runs improve reliability (paper used 64)."
    )
    args = parser.parse_args()

    proteins_file = Path(args.proteins_file)
    if not proteins_file.exists():
        raise FileNotFoundError(f"Proteins file not found: {proteins_file}")

    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = RUNS / proteins_file.stem
    output_dir.mkdir(parents=True, exist_ok=True)

    pickle_path = output_dir / "results.pickle"
    csv_path    = output_dir / "results.csv"

    # Read input protein IDs
    with open(proteins_file, 'r') as f:
        input_proteins_id = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(input_proteins_id)} protein IDs from {proteins_file}")

    proteins = load_protein_database(PROTEIN_DB, input_proteins_id)
    print(f"Running search on {len(proteins)} proteins, k={args.k}, runs={args.runs}")

    pairs = [(k, run) for k in args.k for run in range(args.runs)]

    pool = mp.Pool()
    compute_with_args = partial(compute, proteins=proteins,
                                input_proteins_id=input_proteins_id,
                                center=None)
    results = pool.map(compute_with_args, pairs)

    # Always save pickle
    with open(pickle_path, 'wb') as f:
        pickle.dump(results, f)
    print(f"Pickle results saved to {pickle_path}")

    # Always save CSV
    results_to_csv(results, input_proteins_id, csv_path)


if __name__ == '__main__':
    main()