# # Script to load the Regional Earth System Data Lab
#
#
# #### There are four different versions for the Regional ESDL (RegESDL) at the JupyterHub. They vary in spatial resolution (i.e. 0.0083deg or low 0.083deg) and data storage for temporal or spatial analysis (www.https://www.earthsystemdatalab.net/index.php/interact/data-lab/).
#
# - High and low spatial resolution: Cubes can have different spatial resolution (high 0.0083deg or low 0.083deg). Lower resolution facilitates code development and testing. Higher resolution is recommended for final analysis. The lower resolution RegESDL is a resampled of the higher one.
# - Temporal and spatial data storage: Depening on how data is stored temporal or spatial analysis speed can be improved. We offer both setting, so it is the users choise.
#
#
# Estupinan-Suarez, Lina M.
# 2021.03.04

using ESDL

using ESDLPlots

## This is a high resolution RegESDL for temporal analysis. Observe the "0.0083deg" in the zar file name
chigh = Cube("/home/jovyan/work/datacube/ESDCv2.0.0/esdc-8d-0.0083deg-184x60x60-2.0.0_colombia.zarr/")


## This is a low resolution RegESDL for spatial analysis. Observe the "0.083deg" in the zar file name.
clow = Cube("/home/jovyan/work/datacube/ESDCv2.0.0/esdc-8d-0.083deg-1x336x276-2.0.0_colombia.zarr/")

gpp = chigh[variable="GPP"]
evi = clow[variable="EVIterra"]

plotMAP(gpp, dmin=0, dmax=8, time=Date(2001-01-05))


plotMAP(evi, dmin=0, dmax=1, time=Date(2001-01-05))
