# ## Main figures
# 
# ## Figure 5: Histograms of variance explained by the first three principal components (PCs)¶
# 
# 
# ### Estupinan-Suarez, et al.  (2021). A Regional Earth System Data Lab for Understanding Ecosystem Dynamics: An Example from Tropical South America. Front. Earth Sci. 9:613395. doi: 10.3389/feart.2021.613395
# #### Correspondence to: lestup@bgc-jena.mpg.de, linamaesu@gmail.com
# #### GitHub repository: https://github.com/linamaes/Regional_ESDL
# 
# 
# This script does the following:
# - Plots histograms of PCs explained variance by land cover types
# 
# About the notebook:
# - It is written in Julia 1.3
# - "#" comments in the code are intended to explain specific aspects of the coding
# - New steps in workflows are introduced with bold headers
# 
# March 2021, Max Planck Institute for Biogeochemistry, Jena, Germany

using OnlineStats

using ESDL

using ESDLPlots

using CSV

using Dates

using NPZ

using NetCDF

using OnlineStats

using Plots

using DelimitedFiles

gr(size=(600,400))
default(fmt = :png)

pcasigin = loadCube("pcaComSigLoadings1km2014_qamask")

pcasig = pcasigin[info_PCA="PC_Significance"]

plotMAP(pcasig[Components_or_Var ="sigPC_var1"])

pathin1 = "/my_path_in/folder1/"

pathin2 = "/my_path_in/folder2/"

clcin = ncread(string(pathin1, "/LC_ESA/lc_esa2014.nc"),"Land.cover.class.defined.in.LCCS");

clcmis = map(x->x==210.0 ? missing : x, clcin);

pcasig.axes[1]

clc = CubeMem(CubeAxis[pcasig.axes[1], pcasig.axes[2]], clcmis)

plotMAP(clc)

# Load land cover legend as dictionary
clcdic = include(string(pathin2,"/ESA_LC/legend/ESALClegend2014.jl"));

# Assign dictionary as propertes of the clc cube
clc.properties["labels"]=include(string(pathin2,"/ESA_LC/legend/ESALClegend2014.jl"));

plotMAP(clc)

# Get unique values of land cover types in the study area
clcid = unique(clc[:,:]);

i = 2
clcdic[clcid[i]]

i=2
clcmask = map(x-> ismissing(x) ? missing : x==clcid[i] ? 1 : missing, clc) #[lon=lon,lat=lat]

plotMAP(clcmask)

 pcasig[Components_or_Var="sigPC_var1"]

clcmask

commask1 = map((x,y)->x*y, pcasig[Components_or_Var="sigPC_var1"], clcmask)
commask2 = map((x,y)->x*y, pcasig[Components_or_Var="sigPC_var2"], clcmask)
commask3 = map((x,y)->x*y, pcasig[Components_or_Var="sigPC_var3"], clcmask)

heatmap(commask1[:,:]'[end:-1:1,:], aspect_ratio = :equal)#, dmin=0)

function plotPCshisto(i, clcid, clc, pcasig, titlen; legendYN = :false)
# titlen = clcdic[clcid[i]]

# Generate mask for the selected land cover type
clcmask = map(x-> ismissing(x) ? missing : x==clcid[i] ? 1 : missing, clc)

# Select data of each PCs based on the land cover type (mask generation)
@show(i)
commask1 = map((x,y)->x*y, pcasig[Components_or_Var="sigPC_var1"], clcmask)
commask2 = map((x,y)->x*y, pcasig[Components_or_Var="sigPC_var2"], clcmask)
commask3 = map((x,y)->x*y, pcasig[Components_or_Var="sigPC_var3"], clcmask)

# Filter data that is higher than zero
h1 = filter(x->x>0, skipmissing(commask1[:,:]))
h2 = filter(x->x>0, skipmissing(commask2[:,:]))
h3 = filter(x->x>0, skipmissing(commask3[:,:]));

varname = "Number of pixels"
xname = "Variance explained"

# Plot histograms
hbar = plot(h1, seriestype=:barhist, color="red", linealpha=0.0, label="Component 1",
    title = titlen, titlefont = font(10), legend = legendYN,
        xlabel = xname, xlim = (0,1), xtickfont = font(10),
        ylabel = varname, ytickfont=font(10), guidefont=font(10))
plot!(h2, seriestype=:barhist, color="green", linealpha=0.0, label="Component 2")
plot!(h3, seriestype=:barhist,  color="blue", linealpha=0.0, label="Component 3",
    size=(400,300))

return(hbar)

end

i =2
titlen = clcdic[clcid[i]]

titlen = "Tree cover BrEv-co."

# Replace title with a shorter name
plotPCshisto(i, clcid, clc, pcasig, titlen; legendYN = :true)

i = 15 # Shrubland
titlen = clcdic[clcid[i]]

plotPCshisto(i, clcid, clc, pcasig, titlen; legendYN = :false)

i = 14 # Grassland
titlen = clcdic[clcid[i]]

plotPCshisto(i, clcid, clc, pcasig, titlen; legendYN = :false)

i =12
titlen = clcdic[clcid[i]]

# Replace title with a shorter name
titlen = "Herbaceous cover with trees/shrubs"#

plotPCshisto(i, clcid, clc, pcasig, titlen; legendYN = :false)
