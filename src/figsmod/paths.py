from pathlib import Path

ROOT            = Path(__file__).resolve().parents[2]
DATA            = ROOT / "data"
RAW_DATA        = DATA / "raw"
PROCESSED_DATA  = DATA / "processed"
RUNS            = ROOT / "runs"
RESULTS         = ROOT / "results"
TABLES          = RESULTS / "tables"

PROTEIN_DB      = RAW_DATA / "protein_database_1.pickle"
BIOGRID_NET     = RAW_DATA / "biogrid_net.gpickle"

# Create output directories on import so scripts never raise FileNotFoundError
for _d in [RUNS, RESULTS, TABLES]:
    _d.mkdir(parents=True, exist_ok=True)