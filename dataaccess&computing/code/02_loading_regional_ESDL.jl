# # Data access
# # Script 02.02: Loading the Regional Earth System Data Lab
#
#
# ## Estupinan-Suarez, Gans, et al. (2021). A Regional Earth System Data Lab for Understanding Ecosystem Dynamics: An Example from Tropical South America. Front. Earth Sci. 9:613395. doi: 10.3389/feart.2021.613395
# ### Correspondence to: lestup@bgc-jena.mpg.de, linamaesu@gmail.com
# ### GitHub repository: https://github.com/linamaes/Regional_ESDL


# #### There are four different versions for the Regional ESDL (RegESDL). They vary in (i) spatial resolution (i.e. high 0.0083deg or low 0.083deg) and (ii) data storage sttings for temporal or spatial analysis.
#
# - High and low spatial resolution: Cubes can have different spatial resolutions (high 0.0083deg or low 0.083deg). Lower resolution facilitates code development and testing. Higher resolution is recommended for final analysis. The lower resolution RegESDL is a resampled version of the higher one.
# - Temporal and spatial data storage: Depending on how data is stored temporal or spatial analysis speed can be improved. We offer both settings, the user can choose.
#
#
# June 2021, Max Planck Institute for Biogeochemistry, Jena, Germany


using ESDL # This code requires ESDL 0.9.1

using ESDLPlots


chigh = esdc(region = "Colombia", res="high")

## This is a low resolution RegESDL for spatial analysis. Observe the "0.083deg" in the zar file name.
clow = esdc(region = "Colombia", res = "low") #Does not work at the moment until moved back from glacier

gpp = chigh[variable="GPP"]
evi = clow[variable="EVIterra"]

plotMAP(gpp, dmin=0, dmax=8, time=Date(2001-01-05))

plotMAP(evi, dmin=0, dmax=1, time=Date(2001-01-05))
