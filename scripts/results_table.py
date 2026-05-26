import pickle
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from figsmod.paths import PROTEIN_DB, RUNS, TABLES
from figsmod.protein import Protein


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate per-network pickle files into a single results CSV.",
        usage="python results_table.py [--runs_dir DIR] [--output FILE]",
        epilog="Example: python results_table.py\n"
               "         python results_table.py --runs_dir runs/proteome/ --output results/tables/proteome.csv"
    )
    parser.add_argument("--runs_dir", type=str, default=str(RUNS / "proteome"),
                        help="Directory of per-network pickle files")
    parser.add_argument("--output", type=str,
                        default=str(TABLES / "results_table.csv"),
                        help="Output CSV path")
    parser.add_argument("--no_filter", action="store_true",
                        help="Include unconverged runs and dropped-out proteins")
    args = parser.parse_args()

    runs_dir   = Path(args.runs_dir)
    output     = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    apply_filter = not args.no_filter

    with open(PROTEIN_DB, 'rb') as f:
        protein_database = pickle.load(f)

    data_rows = []
    results = list(runs_dir.glob("*.pickle"))
    print(f"Found {len(results)} result files in {runs_dir}")

    for n, path in enumerate(results, 1):
        print(f"\r{n}/{len(results)}", end='', flush=True)
        protein_hub = path.stem
        with open(path, 'rb') as f:
            data = pickle.load(f)
        for run in data:
            if apply_filter and run['completed_iterations'] == 1000:
                continue
            for i, (protein_id, motif) in enumerate(run['motifs'].items()):
                if apply_filter and run['weights'][i] == 0:
                    continue
                start  = run['positions'][i]
                length = len(motif)
                end    = start + length
                info_content = run['avg_score']

                protein = protein_database.get(protein_id)
                if protein is not None:
                    plm          = float(np.mean(protein.plm[start:end]))
                    disorder     = float(np.mean(protein.disorder[start:end]))
                    solvent_acc  = float(np.mean(protein.rsa[start:end]))
                    conservation = float(np.mean(protein.rlc[start:end]))
                else:
                    plm = disorder = solvent_acc = conservation = float('nan')

                data_rows.append([
                    protein_id, protein_hub, motif,
                    start + 1, end, length, info_content,
                    plm, disorder, solvent_acc, conservation
                ])
    print()

    df = pd.DataFrame(data_rows, columns=[
        'protein', 'binding_protein', 'motif',
        'start', 'end', 'length', 'info_content',
        'plm', 'disorder', 'solvent_acc', 'conservation'
    ])
    df = df.sort_values('protein')

    alpha = 0.6
    w = (5 * df['plm'] + df['disorder'] + df['solvent_acc'] + df['conservation']) / 8
    df['raw_score']        = (alpha * df['info_content'] + (1 - alpha) * w).round(6)
    df['percentile_score'] = df['raw_score'].rank(pct=True).round(6)
    df['info_content']     = df['info_content'].round(4)

    df.to_csv(output, index=False)
    print(f"Results table saved to {output} ({len(df)} rows)")


if __name__ == '__main__':
    main()