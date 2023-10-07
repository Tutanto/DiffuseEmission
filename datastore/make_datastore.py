import os
## plotting libraries

# python library that manages paths
from pathlib import Path

## general mathematical and astronomical libraries
from astropy import units as u
from astropy.coordinates import SkyCoord

## basic imports from Gammapy to manage data 
from gammapy.data import DataStore, EventList, GTI
from gammapy.maps import  MapAxis, WcsGeom
from gammapy.datasets import  MapDataset
from gammapy.makers import MapDatasetMaker, SafeMaskMaker

from GPyUtils.logging_config import logging_conf
from modules.variables import *

# let's read all the event files into the DC folder:
# Name of the output directory where events are
events_dir = dc_folder / "sim_events"

# Read all the new selected event files into the "new_events" folder:
paths = list(events_dir.rglob("events_*.fits"))
paths.sort()


#Here I create a new folder to store the selected events, and select all the events with MC_ID 10001 or 10018
# 10001 --> diffuse, 10018 --> fermi_bubble
# uncomment only one:
IDs, file_name = [10017, 10018], 'fermi'
#IDs, file_name = [10000, 10001], 'diffuse'

path_to_selected_events = path_to_new_events / file_name
(path_to_selected_events).mkdir(parents=True, exist_ok=True)

for i, p in enumerate(paths):
    filenumber = "{0:04d}".format(i)
    events = EventList.read(p)
    gti = GTI.read(p)
    selected_mcid = events.select_parameter(parameter="MC_ID", band=(IDs[0], IDs[1]))
    selected_mcid.write(path_to_selected_events / f"selected_events_{file_name}_{filenumber}.fits", gti=gti, overwrite=True)

# let's read all the new selected event files into the "new_events" folder:
new_paths = list(path_to_selected_events.rglob("selected_events*.fits"))

#let's store the files into a `DataStore` object
data_store = DataStore.from_events_files(events_paths=new_paths, 
                                         irfs_paths=astri_irf)
print(data_store)
# Save the DataStore information to disk
path_to_datastore = path_to_datastores / file_name
(path_to_datastore).mkdir(parents=True, exist_ok=True)
data_store.hdu_table.write(path_to_datastore / f"hdu-index_{file_name}.fits.gz", overwrite=True)
data_store.obs_table.write(path_to_datastore / f"obs-index_{file_name}.fits.gz", overwrite=True)
