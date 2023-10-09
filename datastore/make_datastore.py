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
#IDs, file_name = [10017, 10018], 'fermi'
IDs, file_name = [10001, 10002], 'diffuse'

# Set up log configuration and create a logger
logger = logging_conf(path_to_logs, f"create_{file_name}_datastore.log")

path_to_selected_events = path_to_new_events / file_name
(path_to_selected_events).mkdir(parents=True, exist_ok=True)

logger.debug(f"Saving new events in: {path_to_selected_events}")

for i, p in enumerate(paths):
    filenumber = "{0:04d}".format(i)
    events = EventList.read(p)
    gti = GTI.read(p)
    selected_mcid = events.select_parameter(parameter="MC_ID", band=(IDs[0], IDs[1]))
    selected_mcid.write(path_to_selected_events / f"selected_events_{file_name}_{filenumber}.fits", gti=gti, overwrite=True)
    logger.debug(f"Event saved: {filenumber}")

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
logger.debug(f"Datastore saved in: {path_to_datastore}")

# Load the DataStore from the previously saved directory
data_store = DataStore.from_dir(path_to_datastore, 
                                hdu_table_filename=f"hdu-index_{file_name}.fits.gz", 
                                obs_table_filename=f"obs-index_{file_name}.fits.gz")

# Print information about the DataStore
print(data_store)

# Get a list of observations from the DataStore
table = data_store.obs_table
pos_obs = SkyCoord(table['GLON_PNT'], table['GLAT_PNT'], frame='galactic', unit='deg')
observations = data_store.get_observations()

# Define the path to the dataset directory where analysis results will be saved
path_to_dataset = path_to_datasets / 'width_22x10'
(path_to_dataset).mkdir(parents=True, exist_ok=True)

# Define the name of the dataset file to be created
e_min = 1
e_max = 200
bin = 20
binsz = 0.02
dataset_name = f"dataset_{file_name}_ene_{e_min}_{e_max}_bin_{bin}_binsz_{binsz}.fits.gz"

# Check if the dataset file already exists, and if so, raise an exception
if (path_to_dataset / dataset_name).is_file():
    logger.error(f"{dataset_name} already exists!")
    raise Exception(f"{dataset_name} already exists!")

# Center of the region
pos_target = SkyCoord(73.0, 2.0, frame='galactic', unit='deg')

# Define the energy axis for the dataset
energy_axis = MapAxis.from_energy_bounds(e_min, e_max, bin, unit="TeV")

# Create a WCS geometry for the dataset
geom = WcsGeom.create(
    skydir=pos_target,
    axes=[energy_axis],
    width=[22 * u.deg, 10 * u.deg],
    binsz=binsz * u.deg,
    frame="galactic",
)

# Define the true energy axis for the dataset (used for reduced IRFs)
energy_axis_true = MapAxis.from_energy_bounds(
    0.5, 300, 25, unit="TeV", name="energy_true"
)

# Create an empty stacked MapDataset
stacked = MapDataset.create(
    geom=geom, energy_axis_true=energy_axis_true, name=f"dc_cygnus_{file_name}-stacked"
)

# Create instances of MapDatasetMaker, SafeMaskMaker, and FoVBackgroundMaker
# to process and prepare data for the stacked dataset
maker = MapDatasetMaker(selection=["counts", "background", "exposure", "edisp", "psf"])
maker_safe_mask = SafeMaskMaker(methods=["aeff-default"])

# Log information about the dataset configuration and processing steps
logger.debug(f"Dataset energy: {energy_axis}")
logger.debug(f"Dataset geometry: {geom.wcs}")
logger.debug(f"{maker}")
logger.debug(f"{maker_safe_mask}")

# Loop over each observation, process the data, and stack it onto the final dataset
for obs in observations:
    # First, produce a cutout of the target map
    cutout = stacked.cutout(
        obs.pointing_radec, width="12 deg"
    )
    # Fill the MapDataset with data in this cutout geometry
    dataset = maker.run(cutout, obs)
    # Apply the data quality cut using SafeMaskMaker
    dataset = maker_safe_mask.run(dataset, obs)
    # Stack the resulting dataset cutout onto the final dataset
    stacked.stack(dataset)

# Log the final stacked dataset information
logger.info(stacked)

# Save the stacked DataSet to disk
filename = path_to_dataset / dataset_name
stacked.write(filename, overwrite=False)
logger.debug(f"Dataset saved: {filename}")
logger.debug("#### END ####")
