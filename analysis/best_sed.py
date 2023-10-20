# Import necessary libraries and modules
import numpy as np
import astropy.units as u

from gammapy.datasets import MapDataset
from gammapy.modeling.models import Models
from gammapy.estimators import FluxPointsEstimator

# Import custom modules and functions from GPyUtils and other sources
from GPyUtils.general_utils import calculate_sed
from GPyUtils.logging_config import logging_conf
from modules.variables import *

# Name of the dataset file to be used for the SED calculation
diffuse = 'no_diffuse'
file_name = 'all_IDs'
strategy = 1
tol = 0.01
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
path_to_fit_models = path_to_results / 'spatial_freeze' / f'strategy_{strategy}' / f'tol_{tol}' / diffuse / f"{file_name}_ene_{e_min}_{e_max}_bin_{bin}_binsz_{binsz}" / "models"
path = list(path_to_fit_models.rglob("01*.yaml"))
# Define the names of the output files for the fit results
filename = f"{path[0].with_suffix('').name}_fitted"

# Load the fitted model from the YAML file
models_fit = Models.read(path[0])

# Set the loaded model as the source model in the dataset
dataset.models = models_fit

# Get the first model in the dataset, assuming only one source is present
source = dataset.models[0]

# Define the name of the output directory for the SED data points
path_to_result_fluxpoints = path_to_fluxdatapoints / 'spatial_freeze' / f'strategy_{strategy}' / f'tol_{tol}' / diffuse / f"{file_name}_ene_{e_min}_{e_max}_bin_{bin}_binsz_{binsz}"
path_to_result_fluxpoints.mkdir(parents=True, exist_ok=True)
saved_fluxpoints = path_to_result_fluxpoints / f'{filename}.fits'

# Set up log configuration and create a logger for the SED calculation
logger = logging_conf(path_to_logs, f"best_sed_{filename}.log")

# Log information about the dataset and model used for the SED calculation
logger.debug(f"Spectrum extraction using dataset: {dataset.name}")
logger.debug(f"Dataset geom: {dataset.geoms['geom']}")
logger.debug(f"Dataset energy axis: {dataset.geoms['geom'].axes['energy']}")

# Compute fluxpoints
sed_points = np.logspace(0, 2, num=10) * u.TeV
fpe = FluxPointsEstimator(energy_edges=sed_points, source="cygnus_diffuse", selection_optional=["ul"])
# Compute the flux data points for the specified source and energy range
flux_points = fpe.run(datasets=[dataset])
flux_points.write(filename=saved_fluxpoints, overwrite=True)

#result = calculate_sed(dataset=dataset, source=source, emin=1, emax=100, step=2., n_jobs=20, logger=logger)
# Save the calculated SED data points to a FITS file
#result.write(filename=output_flux / f"{sedname}.fits", overwrite=True)

# Log information about the location of the saved SED file
logger.info(f"Spectrum saved in: {saved_fluxpoints}")
