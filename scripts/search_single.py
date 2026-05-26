import pickle
import argparse
import multiprocessing as mp
from functools import partial
from pathlib import Path
from figsmod.paths import PROTEIN_DB, BIOGRID_NET, RUNS
from figsmod.search import read_biogrid, load_protein_database, compute
from figsmod.protein import Protein


def main():
    parser = argparse.ArgumentParser(
        description="Run FiGS-MoD on a single LMBD protein network from BioGRID.",
        usage="python search_single.py <lmbd_protein> [output_path]",
        epilog="Example: python search_single.py Q15084\n"
               "         python search_single.py Q15084 runs/Q15084.pickle"
    )
    parser.add_argument("center", type=str,
                        help="UniProt ID of the LMBD hub protein")
    parser.add_argument("output_path", type=str, nargs="?", default=None,
                        help="Output pickle path (default: runs/<center>.pickle)")
    args = parser.parse_args()

    output_path = Path(args.output_path) if args.output_path else RUNS / f"{args.center}.pickle"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    num_run = 1
    ks = [7]

    input_proteins_id = read_biogrid(BIOGRID_NET, args.center)
    proteins = load_protein_database(PROTEIN_DB, input_proteins_id)
    pairs = [(k, run) for k in ks for run in range(num_run)]

    pool = mp.Pool()
    compute_with_args = partial(compute, proteins=proteins,
                                input_proteins_id=input_proteins_id,
                                center=args.center)
    results = pool.map(compute_with_args, pairs)

    with open(output_path, 'wb') as f:
        pickle.dump(results, f)
    print(f"Results saved to {output_path}")


if __name__ == '__main__':
    main()