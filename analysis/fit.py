# Import necessary libraries and modules
import numpy as np
import astropy.units as u

from gammapy.modeling import Fit
from gammapy.datasets import MapDataset
from gammapy.estimators import FluxPointsEstimator
from gammapy.modeling.models import FoVBackgroundModel, Models

# Import custom modules and functions from GPyUtils and other sources
from GPyUtils.logging_config import logging_conf
from GPyUtils.write_minuit_result import MinuitResultWriter
from modules.variables import *

# Define the name of the dataset file to be used for the fit
diffuse = 'diffuse'
file_name = 'all_IDs'
strategy = 1
tol = 0.01
e_min = 1
e_max = 200
bin = 20
binsz = 0.02
dataset_name = f"dataset_{file_name}_ene_{e_min}_{e_max}_bin_{bin}_binsz_{binsz}.fits.gz"

# Set the path to the dataset file
path_to_dataset = path_to_datasets / 'width_22x10' / dataset_name

# Reading the Dataset from the specified file
dataset = MapDataset.read(filename=path_to_dataset)

# Define the filename and path for the fitted model to be saved
models = path_to_models / diffuse
saved_models = path_to_results / 'multi_models' / f'strategy_{strategy}' / f'tol_{tol}' / diffuse / f"{file_name}_ene_{e_min}_{e_max}_bin_{bin}_binsz_{binsz}" / "models"
saved_jsons = path_to_results / 'multi_models' / f'strategy_{strategy}' / f'tol_{tol}' / diffuse / f"{file_name}_ene_{e_min}_{e_max}_bin_{bin}_binsz_{binsz}" / "jsons"

saved_models.mkdir(parents=True, exist_ok=True)
saved_jsons.mkdir(parents=True, exist_ok=True)

paths = list(models.rglob("*.yaml"))
paths.sort()

for path in paths:
    try:
        models_fit = Models.read(path)
        bkg_model = FoVBackgroundModel(dataset_name=dataset.name)
        models_fit.insert(-1, bkg_model)

        dataset.models = models_fit
        
        # Define the name of the source
        if path.with_suffix('').name == 'nullhypothesis':
            model_name = None
        elif path.with_suffix('').name == '00_template':
            model_name = ['fermi bubble']
        else:
            model_name = ['cygnus_diffuse']

        # Define the names of the output files for the fit results
        filename = f"{path.with_suffix('').name}_fitted"
        result_model = f'{filename}.yaml'
        result_json = f'{filename}.json'

        # Set up log configuration and create a logger for the fit
        logger = logging_conf(path_to_logs, f"fit_{filename}.log")

        # Check if the output files already exist, and if so, raise an exception
        if (saved_models / result_model).is_file() or (saved_jsons / result_json).is_file():
            logger.error(f"{result_model} or {result_json} already exists!")
            raise Exception(f"{result_model} or {result_json} already exists!")
        else:
            # Log information about the dataset and model used for the fit
            logger.debug(f"Fitting using dataset: {dataset_name}")
            logger.debug(f"Dataset geom: {dataset.geoms['geom']}")
            logger.debug(f"Dataset energy axis: {dataset.geoms['geom'].axes['energy']}")
            logger.debug(f"Model used: {result_model}")

            # Run the fit using the Minuit optimizer from Gammapy
            fit = Fit(store_trace=True)
            minuit_opts = {"tol": tol, "strategy": strategy}
            fit.backend = "minuit"
            fit.optimize_opts = minuit_opts
            result_fit = fit.run(dataset)

            # Log the fit results
            logger.info(result_fit)

            # Save the fitted model to a YAML file
            dataset.models.write(saved_models / result_model, overwrite=True)
            logger.info(f"Fitted Model saved: {saved_models / result_model}")

            # If the fit was done using the Minuit backend, save the Minuit results to a JSON file
            if fit.backend == "minuit":
                # Create an instance of MinuitResultWriter to handle the results
                stored_results = MinuitResultWriter(
                    dataset,
                    dataset_name,
                    filename,
                    model_name,
                    result_fit
                )

                # Save the Minuit results to a JSON file
                stored_results.save_result_json(saved_jsons / result_json)
                logger.info(f"Minuit results saved: {saved_jsons / result_json}")

                # to print the Minuit results in a readable format.
                # stored_results.pretty_print_result()
            
            if path.with_suffix('').name != 'nullhypothesis':
                # Compute fluxpoints
                sed_points = np.logspace(0, 2, num=11) * u.TeV
                fpe = FluxPointsEstimator(energy_edges=sed_points, source=model_name[0], selection_optional=["ul"])
                # Log information about the dataset and model used for the SED calculation
                logger.debug(f"Spectrum extraction using dataset: {dataset.name}")
                logger.debug(f"Dataset geom: {dataset.geoms['geom']}")
                logger.debug(f"Dataset energy axis: {dataset.geoms['geom'].axes['energy']}")
                logger.debug(f"Computing SED from {sed_points.min} to {sed_points.max()} TeV, with {len(sed_points)} bins")

                flux_points = fpe.run(datasets=[dataset])
            
                # Save fluxpoints
                path_to_result_fluxpoints = path_to_fluxdatapoints / 'multi_models' / f'strategy_{strategy}' / f'tol_{tol}' / diffuse / f"{file_name}_ene_{e_min}_{e_max}_bin_{bin}_binsz_{binsz}"
                path_to_result_fluxpoints.mkdir(parents=True, exist_ok=True)
                saved_fluxpoints = path_to_result_fluxpoints / f'{filename}.fits'
                flux_points.write(filename=saved_fluxpoints, overwrite=True)
                # Log information about the location of the saved SED file
                logger.info(f"Spectrum saved in: {saved_fluxpoints}")
    
    except Exception as e:
        # Log the exception, print an error message, and continue to the next iteration
        logger.error(f"An error occurred in the loop iteration for {path}: {e}")
        continue
