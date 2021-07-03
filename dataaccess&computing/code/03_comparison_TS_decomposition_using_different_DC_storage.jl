# ## Computational Performance
#
# ## Script 03: Computational comparison of time series decomposition using different Data Cubes storage
#
# ### Estupinan-Suarez, Gans, et al. (2021). A Regional Earth System Data Lab for Understanding Ecosystem Dynamics: An Example from Tropical South America. Front. Earth Sci. 9:613395. doi: 10.3389/feart.2021.613395
# #### Correspondence to: lestup@bgc-jena.mpg.de, linamaesu@gmail.com
# #### GitHub repository: https://github.com/linamaes/Regional_ESDL
#
#
# This script illustrates the time difference on time series computing when using different cube storage versions. We provide two examples:
#
# - Example 1: Time series are decomposed in different time-scales i.e. Longer-term, seasonal cycle and short(fast) oscillation, plus the trend
# - Example 2: It finds the maximum value along a time period at pixel level
#
#
#
# About the notebook
# - It is written in Julia 1.3
# - "#" comments in the code are intended to explain specific aspects of the coding
# - New steps in the workflow are introduced with bold headers
#
# March 2021, Max Planck Insitute for Biogeochemistry, Jena, Germany

using ESDL

using ESDLPlots

pathin = "../mypath/"# => Set your path

# Cube storage version for temporal analysis
ctem = Cube(string(pathin,"/Cube_2020highColombiaCube_184x120x120.zarr/"))

# Cube storage version for spatial analysis
cspa = Cube(string(pathin,"/Cube_2020highColombiaCube_1x3360x2760.zarr/"))

yrstart = 2001
yrend = 2005
latsub =(0,1)

@time ctemfft = filterTSFFT(ctem[time=yrstart:yrend, lat=latsub])

@time cspafft = filterTSFFT(cspa[time=yrstart:yrend, lat=latsub])

@time ctempax= mapslices(maximum, ctem[time=yrstart:yrend, lat=latsub], dims="time")

@time cspamax = mapslices(maximum, cspa[time=yrstart:yrend, lat=latsub], dims="time")
