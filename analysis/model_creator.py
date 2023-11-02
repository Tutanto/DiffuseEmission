import itertools
import astropy.units as u

from gammapy.modeling.models import (
    ExpCutoffPowerLawSpectralModel,
    BrokenPowerLawSpectralModel,
    LogParabolaSpectralModel,
    PowerLawSpectralModel,
    GaussianSpatialModel,
    DiskSpatialModel,
    SkyModel,
    Models
)
from modules.variables import *

diffuse = False

original_model_path = dc_folder / 'skymodel' / 'global_skymodel_v1.0.yaml'
original_model = Models.read(original_model_path)

original_model.freeze('spatial')
original_model.freeze('spectral')

for model in original_model:
        if len(model.spectral_model.parameters.names) > 1 and model.spectral_model.tag[0] != 'PowerLawNormSpectralModel':
            model.spectral_model.amplitude.frozen = False
        else:
            model.spectral_model.norm.frozen = False


if not diffuse:
    original_model.remove('diffuse')
    folder_name = 'no_diffuse'
elif diffuse:
    folder_name = 'diffuse'


path_to_original_model = path_to_models / folder_name / '00_template.yaml'
original_model.write(path_to_original_model, overwrite=True)
original_model.remove('fermi bubble')


path_to_original_model = path_to_models / folder_name / 'nullhypothesis.yaml'
original_model.write(path_to_original_model, overwrite=True)

disk = DiskSpatialModel(
    lon_0=79.6188337 * u.deg,
    lat_0=0.962326 * u.deg,
    r_0="3 deg",
    e=0.01,
    phi=0 * u.deg,
    edge_width=0.1,
    frame="galactic",
)
gaussian = GaussianSpatialModel(
    lon_0=79.6188337 * u.deg,
    lat_0=0.962326 * u.deg,
    sigma="2.1 deg",
    e=0.01,
    phi=0 * u.deg,
    frame="galactic",
)
pwl = PowerLawSpectralModel(
        index=2.5,
        amplitude="1e-11 TeV-1 cm-2 s-1",
        reference=1 * u.TeV,
    )
brk = BrokenPowerLawSpectralModel(
    index1=2.5,
    index2=2.7,
    amplitude="1e-13 TeV-1 cm-2 s-1",
    ebreak="10 TeV",
)
log = LogParabolaSpectralModel(
    alpha=2.5,
    amplitude="4e-11 cm-2 s-1 TeV-1",
    reference=1 * u.TeV,
    beta=0.07,
)
cut = ExpCutoffPowerLawSpectralModel(
    amplitude=3e-11 * u.Unit("cm-2 s-1 TeV-1"),
    index=2.5,
    lambda_=0.01 * u.Unit("TeV-1"),
    reference=1 * u.TeV,
)

center_position = ['center_fixed', 'center_free']
spectral_models = [pwl, cut, log, brk]
spatial_models = [disk, gaussian]

for i, (position, spectral, spatial) in enumerate(itertools.product(center_position, spectral_models, spatial_models)):
    j = i+1
    diffuse = SkyModel(spectral_model=spectral, spatial_model=spatial, name="cygnus_diffuse")

    models_diffuse = original_model.copy()
    models_diffuse.insert(0, diffuse)
    if position == "center_fixed":
        models_diffuse[0].spatial_model.lon_0.frozen = True
        models_diffuse[0].spatial_model.lat_0.frozen = True
    if models_diffuse[0].spatial_model.tag[0] == disk.tag[0]:
        models_diffuse[0].parameters["r_0"].max = 4
        models_diffuse[0].parameters["r_0"].min = 0.05
    if models_diffuse[0].spatial_model.tag[0] == gaussian.tag[0]:
        models_diffuse[0].parameters["sigma"].max = 3
        models_diffuse[0].parameters["sigma"].min = 0.1        
    if models_diffuse[0].spectral_model.tag[0] == brk.tag[0]:
        models_diffuse[0].parameters["index1"].max = 3
        models_diffuse[0].parameters["index1"].min = 2
        models_diffuse[0].parameters["index2"].max = 3
        models_diffuse[0].parameters["index2"].min = 2.1
    if models_diffuse[0].spectral_model.tag[0] == cut.tag[0]:
        models_diffuse[0].parameters["lambda_"].max = 0.1
        models_diffuse[0].parameters["lambda_"].min = 0.001
        models_diffuse[0].parameters["index"].max = 3
        models_diffuse[0].parameters["index"].min = 2.1
    if models_diffuse[0].spectral_model.tag[0] == log.tag[0]:
        models_diffuse[0].parameters["alpha"].max = 3
        models_diffuse[0].parameters["alpha"].min = 2
        models_diffuse[0].parameters["beta"].max = 1
        models_diffuse[0].parameters["beta"].min = 0.01
        
    print(models_diffuse)

    trimmed_string_spec = models_diffuse[0].spectral_model.tag[1]
    trimmed_string_spat = models_diffuse[0].spatial_model.tag[1]

    path_to_model = path_to_models / folder_name / f'{j:02d}_{trimmed_string_spec}_{trimmed_string_spat}_{position}.yaml'
    models_diffuse.write(path_to_model, overwrite=True)
