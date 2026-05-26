import subprocess
import pickle
import argparse
from pathlib import Path
from figsmod.paths import BIOGRID_NET, RUNS


def main():
    parser = argparse.ArgumentParser(
        description="Run FiGS-MoD on all LMBD protein networks in BioGRID.",
        usage="python run_proteome.py [output_folder]",
        epilog="Example: python run_proteome.py runs/proteome_20250525/"
    )
    parser.add_argument("output_folder", type=str, nargs="?",
                        default=str(RUNS / "proteome"),
                        help="Directory to save per-network pickle files")
    args = parser.parse_args()

    results_folder = Path(args.output_folder)
    results_folder.mkdir(parents=True, exist_ok=True)

    with open(BIOGRID_NET, 'rb') as f:
        biogrid = pickle.load(f)

    centers = []
    for center in biogrid.nodes:
        if center == '-':
            continue
        filtered = set()
        if len(biogrid[center]) <= 50 and len(biogrid[center]) >= 10:
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
        centers.append(center)  # was centers.append[center] — TypeError fixed

    search_script = Path(__file__).parent / "search_single.py"
    for center in centers:
        out = results_folder / f"{center}.pickle"
        try:
            subprocess.run(
                ['python', str(search_script), center, str(out)],
                check=True
            )
        except Exception as e:
            print(f"Error at {center}: {e}")


if __name__ == '__main__':
    main()