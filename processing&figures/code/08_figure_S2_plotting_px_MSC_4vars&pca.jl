# ## Supplementary figures
# ## Figure S2: Mean Seasonal Cycle (MSC) comparison of vegetation variables and the first principal component (PC) by pixels
# 
# 
# ### Estupinan-Suarez, et al. (2021). A Regional Earth System Data Lab for Understanding Ecosystem Dynamics: An Example from Tropical South America. Front. Earth Sci. 9:613395. doi: 10.3389/feart.2021.613395
# #### Correspondence to: lestup@bgc-jena.mpg.de, linamaesu@gmail.com
# #### GitHub repository: https://github.com/linamaes/Regional_ESDL
# 
# This script does the following:
# Steps:
# - Loads vegetation variables and calculates their MSC
# - Loads computed MSC from the first PCA component
# - Standardize data for plotting purpose
# - Plots vegetation variables MSC at pixel level for different land cover types 
# 
# 
# About the notebook
# - It is written in Julia 1.3
# - "#" comments in the code are intended to explain specific aspects of the coding
# - New steps in workflows are introduced with bold headers
# 
# April 2021, Max Planck Institute for Biogeochemistry, Jena, Germany

using OnlineStats

using ESDL

using CSV

using Dates

using NPZ

using NetCDF

using Plots

using ESDLPlots

using DelimitedFiles

using MultivariateStats

gr(size=(600,400))
default(fmt = :png)

pathin1 = "/my_pathin1/.../"
pathin2 = "/my_pathin2/.../"
pathout = "/my_pathout/.../"

c = Cube(string(pathin1, "/Cube_2020highColombiaCube_184x120x120.zarr/"))

# set missings
cmis = map(x-> ismissing(x) ? missing : x==255.0 ? missing : x==65335.0 ? missing : x, c)

# select assessed time period
csub = cmis[time=2001:2014, Variable=["gross_primary_productivity", "ndvi",
        "evi","fpar"]]

cpca = loadCube("pcacom1msc1km2014_mask")

clcin = ncread(string(pathin2, "data/LC_ESA/lc_esa2014.nc"),"Land.cover.class.defined.in.LCCS")
clcmis = map(x->x==210.0 ? missing : x, clcin)
clc = CubeMem(CubeAxis[cpca.axes[1], cpca.axes[2]], clcmis)
clcdic = include("lc_legend_abr2.jl");

# standardized function
function stdFx(xout, xin)
    xout[:] .= map(x -> (x - OnlineStats.mean(xin))/sqrt(var(xin)), xin)
end
function stdFx(xout, xin)
    xout[:] .= map(x -> (x - mean(xin))/sqrt(var(xin)), xin)
end
# smoothing function

function getIndex_circ2(x::Vector,i::Integer)
while i<1 || i>length(x)
i = i<1 ? i+length(x) : i-length(x)
end
x[i]
end

function smoothMSC(xin;ns=2)
    any(ismissing,xin) && return(fill!(similar(xin),missing))
    scur = sum(getIndex_circ2(xin,i) for i=1-ns:1+ns)
    afac = 1/(2ns+1)
    xout = zeros(size(xin))
    for i=1:length(xin)
        xout[i]=scur*afac
        scur+=getIndex_circ2(xin,i+ns+1)-getIndex_circ2(xin,i-ns)
    end
    xout
end

# function to Normalize data
function norFx(xout, xin)
    xout[:] .= map(x -> (x - minimum(xin))/(maximum(xin) - ninimum(xin)), xin)
end

# extract longitude and latitude coordinates
lonaxc = collect(csub.cubeaxes[1].values)
lataxc = collect(csub.cubeaxes[2].values);

# create an empty coordinates matrix
pxcor = Array{Float32,2}(undef, 5, 4);

# define pixels location and coordinates
pxcor[:,1] .= [1500, 1830, 1350, 1830, 1801]
pxcor[:,2] .= [1000, 1350, 1830, 2000, 1081]
pxcor[:,3] =  map(x->lonaxc[Int.(x)], (pxcor[:,1]))
pxcor[:,4] =  map(x->lataxc[Int.(x)], (pxcor[:,2]));

# start settings
i = 3
nsx = 5
cin = csub

# get one pixel cube, gap filles and decomposed time series
cpx = cin[lon=pxcor[i,3], lat=pxcor[i,4]]
cfill = gapFillMSC(cpx)
cfft = filterTSFFT(cfill)

# smooth PCA MSC
pxpca = cpca[lon=pxcor[i,3], lat=pxcor[i,4]]
spxpca = smoothMSC(pxpca[:], ns=nsx);

# fing pixel's land cover
iddic = clc[Int(pxcor[i,1]), Int(pxcor[i,2])]
titlex = clcdic[iddic]

# standardized variables
indims = InDims("Time")
outdims = OutDims("Time")
cstd = mapCube(stdFx, cfill, indims = indims, outdims = outdims)

# get MSC and smooth it
cmsc = getMSC(cstd)
csmo = mapslices(smoothMSC, cmsc, dims="MSC", ns=6)

# get days of year from MSC cube axis
taxin = collect(cmsc.axes[1].values)
tax = map(x->Dates.dayofyear(x), taxin[1:46]);

# get coordinates for title
coordout = string("(Lat. ", round(pxcor[i,4], digits=2) ,"°, Lon. ", round(pxcor[i,3], digits=2),"°)")#,

# plot variables
pout = plot(tax, csmo[:,1], label="GPP", lw=2.5,  xlabel="Time (day of year)",  ylabel="Standardized values
of vegetation variables",
         title = string(titlex, "\n", coordout), titlefontsize = 10, legendfontsize = 7,
        legend =  legend=:bottomleft, size=(400,300), colour=:green)
    plot!(tax, csmo[:,2], lw=2, label="NDVI", colour=:purple)
    plot!(tax, csmo[:,3], lw=2, label="EVI", colour=:red)
    plot!(tax, csmo[:,4], lw=2, label="FPAR", colour=:black)
    plot!(tax, spxpca[:]*-1, lw=2, label="PC1", colour=:blue)

# function for plotting vegetation variables of a single pixel

function pxMSCvarsmo(cin, cpca, clc, clcdic, pxcor, i, smoothMSC; nsx=3)

    cpx = cin[lon=pxcor[i,3], lat=pxcor[i,4]]
    cfill = gapFillMSC(cpx)

    pxpca = cpca[lon=pxcor[i,3], lat=pxcor[i,4]]
    spxpca = smoothMSC(pxpca[:], ns=nsx)

    iddic = clc[Int(pxcor[i,1]), Int(pxcor[i,2])]
    titlex = clcdic[iddic]

    indims = InDims("Time")
    outdims = OutDims("Time")
    cstd = mapCube(stdFx, cfill, indims = indims, outdims = outdims)
    cmsc = getMSC(cstd)
    csmo = mapslices(smoothMSC, cmsc, dims="MSC", ns=6)

    coordout = string("(Lat. ", round(pxcor[i,4], digits=2) ,"°, Lon. ", round(pxcor[i,3], digits=2),"°)")#,

    pout = plot(tax, csmo[:,1], label="GPP", lw=2.5,  xlabel="Time (day of year)",  ylabel="Standardized values
of vegetation variables",
         title = string(titlex, "\n", coordout), titlefontsize = 10, legendfontsize = 7,
        legend =  legend=:bottomleft, size=(400,300), colour=:green)
    plot!(tax, csmo[:,2], lw=2, label="NDVI", colour=:purple)
    plot!(tax, csmo[:,3], lw=2, label="EVI", colour=:red)
    plot!(tax, csmo[:,4], lw=2, label="FPAR", colour=:black)
    plot!(tax, spxpca[:]*-1, lw=2, label="PC1", colour=:blue)

    return pout

end

px1 = pxMSCvarsmo(csub, cpca, clc, clcdic, pxcor, 1, smoothMSC)

savefig(px1, string(pathout, "plots/mscpxvarpca/1km/mask/MSCvar&pca2014flooded_mask.png"))

px2 = pxMSCvarsmo(csub, cpca, clc, clcdic, pxcor, 2, smoothMSC, nsx=2)

savefig(px2, string(pathout, "plots/mscpxvarpca/1km/mask/MSCvar&pca2014trees_mask.png"))

px3 = pxMSCvarsmo(csub, cpca, clc, clcdic, pxcor, 3, smoothMSC, nsx=2)

savefig(px3, string(pathout, "plots/mscpxvarpca/1km/mask/MSCvar&pca2014trees2_mask.png"))

px4 = pxMSCvarsmo(csub, cpca, clc, clcdic, pxcor, 4, smoothMSC)

savefig(px4, string(pathout, "plots/mscpxvarpca/1km/mask/MSCvar&pca2014trees3_mask.png"))

px5 = pxMSCvarsmo(csub, cpca, clc, clcdic, pxcor, 5, smoothMSC)

savefig(px5, string(pathout, "plots/mscpxvarpca/1km/mask/MSCvar&pca2014grassland_mask.png"))
