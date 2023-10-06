import os
from pathlib import Path

cwd = Path.cwd()
home = Path.home()
root_dir = cwd.parents[0]

dc_data = Path(os.environ["ASTRI_DATA"], "sim_events")
astri_irf = Path(os.environ["ASTRI_IRF"],"astri_100_43_008_0502_C0_20_AVERAGE_50h_SC_v1.0.lv3.fits")
