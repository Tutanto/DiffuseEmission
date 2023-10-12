# Import necessary libraries and modules
from gammapy.datasets import MapDataset
from gammapy.modeling.models import Models

# Import custom modules and functions from GPyUtils and other sources
from GPyUtils.general_utils import calculate_sed
from GPyUtils.logging_config import logging_conf
from modules.variables import *

# Name of the dataset file to be used for the SED calculation
file_name = 'all_IDs'
e_min = 0.7
e_max = 100
bin = 20
binsz = 0.02
dataset_name = f"dataset_{file_name}_ene_{e_min}_{e_max}_bin_{bin}_binsz_{binsz}.fits.gz"

# Set the path to the dataset file
path_to_dataset = path_to_datasets / 'width_22x10' / dataset_name

# Reading the MapDataset from the specified file
dataset = MapDataset.read(filename=path_to_dataset)

# Name of the YAML file containing the fitted model to be used for the SED calculation
filename = "00_pl_disk_center_fixed_fitted"
path_to_fit_models = path_to_results / 'no_diffuse' / f"{file_name}_ene_{e_min}_{e_max}_bin_{bin}_binsz_{binsz}" / "models"

# Load the fitted model from the YAML file
models_fit = Models.read(f"{path_to_fit_models}/{filename}.yaml")

# Set the loaded model as the source model in the dataset
dataset.models = models_fit

# Get the first model in the dataset, assuming only one source is present
source = dataset.models[0]

# Define the name of the output directory for the SED data points
sedname = f"spectrum_{source.name}"
output_flux = path_to_fluxdatapoints / filename
output_flux.mkdir(parents=True, exist_ok=True)

# Set up log configuration and create a logger for the SED calculation
logger = logging_conf(path_to_logs, f"best_sed_{filename}.log")

# Log information about the dataset and model used for the SED calculation
logger.debug(f"Spectrum extraction using dataset: {dataset.name}")
logger.debug(f"Dataset geom: {dataset.geoms['geom']}")
logger.debug(f"Dataset energy axis: {dataset.geoms['geom'].axes['energy']}")

# Compute the flux data points for the specified source and energy range
result = calculate_sed(dataset=dataset, source=source, emin=1, emax=100, step=2., n_jobs=20, logger=logger)

# Save the calculated SED data points to a FITS file
result.write(filename=output_flux / f"{sedname}.fits", overwrite=True)

# Log information about the location of the saved SED file
logger.info(f"Spectrum saved in: {output_flux}")
