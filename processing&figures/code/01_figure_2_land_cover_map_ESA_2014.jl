# ## Main figures
# ## Figure 2: Land cover map for the regional ESDL
#
# ### Estupinan-Suarez, et al. (2021). A Regional Earth System Data Lab for Understanding Ecosystem Dynamics: An Example from Tropical South America. Front. Earth Sci. 9:613395. doi: 10.3389/feart.2021.613395
# #### Correspondence to: lestup@bgc-jena.mpg.de, linamaesu@gmail.com
# #### GitHub repository: https://github.com/linamaes/Regional_ESDL
#
# This script does the following:
# - Loads land cover data from the European Space Agency (ESA) map in 2014
#
# About the notebook
# - It is written in Julia 1.3
# - "#" comments in the code are intended to explain specific aspects of the coding
# - New steps in workflows are introduced with bold headers
#
# March 2021, Max Planck Institute for Biogeochemistry, Jena, Germany

using ESDL

using ESDLPlots

pathin = "/my_path/.../"

clc = loadCube(string(pathin, "ESA_LC/ColombiaLCCube_high.zarr/"))

plotMAP(clc, year=2014)

function unit8dataFx(xout, xin)
    xout .= Int.(xin)
    #return xout
end

indims = InDims("Lon","Lat")
outdims = OutDims("Lon","Lat",outtype=Int)

clcyear = mapCube(unit8dataFx, clc[Year=2014], indims = indims, outdims = outdims)

# Add legend to data
clcyear.properties["labels"]=include(string(pathin,"ESA_LC/legend/ESALClegend.jl"))

# colors differ from the map shown in the refered publication (figure 2)
plotMAP(clcyear)
