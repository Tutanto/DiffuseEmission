import os
from pathlib import Path

cwd = Path.cwd()
home = Path.home()
root_dir = cwd.parents[0]

path_to_logs = root_dir / "logs"
path_to_datasets = root_dir / "datasets"
path_to_datastores = root_dir / "datastore"
path_to_new_events = root_dir / "events"

dc_folder = Path(os.environ["ASTRI_DATA"], "DC-GCygni")
astri_irf = Path(os.environ["ASTRI_IRF"],"astri_100_43_008_0502_C0_20_AVERAGE_50h_SC_v1.0.lv3.fits")
