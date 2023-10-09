from gammapy.data import EventList, GTI, DataStore
from astropy import units as u
from astropy.coordinates import SkyCoord
from gammapy.maps import  MapAxis, WcsGeom
from gammapy.datasets import  MapDataset
from gammapy.makers import MapDatasetMaker, SafeMaskMaker

# Import custom modules and functions from GPyUtils and other sources
from GPyUtils.logging_config import logging_conf
from modules.variables import *

# Set up log configuration and create a logger
logger = logging_conf(path_to_logs, f"combine_datastore.log")

# let's read all the event files:
path_to_diffuse = path_to_new_events / 'diffuse'
path_to_fermi = path_to_new_events / 'fermi'
paths_diff = list(path_to_diffuse.rglob("selected_events*.fits"))
paths_fermi = list(path_to_fermi.rglob("selected_events*.fits"))
logger.debug(f'Read files in {path_to_diffuse}')
logger.debug(f'Read files in {path_to_fermi}')
new_path = path_to_new_events / 'diffuse+fermi'
logger.debug(f"Created {new_path}")
#Here I create a new folder to store the selected events, and select all the events with MC_ID 10001
(new_path).mkdir(parents=True, exist_ok=True)


for i, (filename_1, filename_2) in enumerate(zip(paths_diff, paths_fermi)):
    filenumber = "{0:04d}".format(i)
    events_1, events_2 = EventList.read(filename_1), EventList.read(filename_2)
    gti_1, gti_2 = GTI.read(filename_1), GTI.read(filename_2)
    # stack in place, now the _1 object contains the information of both
    gti_1.stack(gti_2)
    events_1.stack(events_2)
    events_1.write(new_path / f"selected_events_diffuse+fermi_{filenumber}.fits", gti=gti_1, overwrite=True)

# let's read all the new selected event files into the "new_events" folder:
new_paths = list(new_path.rglob("selected_events*.fits"))
#let's store the files into a `DataStore` object
data_store = DataStore.from_events_files(events_paths=new_paths, 
                                         irfs_paths=astri_irf)
print(data_store)
# Let's save them on disk:
data_store.hdu_table.write(path_to_datastores / "selected_diffuse+fermi_hdu-index.fits.gz", overwrite=True)
data_store.obs_table.write(path_to_datastores / "selected_diffuse+fermi_obs-index.fits.gz", overwrite=True)
# Let's read the DataStore
data_store = DataStore.from_dir(path_to_datastores, 
                                "selected_diffuse+fermi_hdu-index.fits.gz", 
                                "selected_diffuse+fermi_obs-index.fits.gz")
data_store.info()
observations = data_store.get_observations()
table = data_store.obs_table
pos_obs = SkyCoord(table['GLON_PNT'], table['GLAT_PNT'], frame='galactic', unit='deg')

# Define the path to the dataset directory where analysis results will be saved
path_to_dataset = path_to_datasets / 'width_22x10' / 'diffuse+fermi'
(path_to_dataset).mkdir(parents=True, exist_ok=True)

# Define the name of the dataset file to be created
e_min = 1
e_max = 200
bin = 20
binsz = 0.02
dataset_name = f"dataset_diffuse+fermi_ene_{e_min}_{e_max}_bin_{bin}_binsz_{binsz}.fits.gz"

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
    geom=geom, energy_axis_true=energy_axis_true, name=f"dc_cygnus_diffuse+fermi-stacked"
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