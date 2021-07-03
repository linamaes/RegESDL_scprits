# ## Main figures
# 
# ## Figure 7: Computation of the mean the seasonal cycle (MSC) and standard deviation (SDSC) by biotic units (bu)
# 
# ### Estupinan-Suarez, et al. (2021). A Regional Earth System Data Lab for Understanding Ecosystem Dynamics: An Example from Tropical South America. Front. Earth Sci. 9:613395. doi: 10.3389/feart.2021.613395
# #### Correspondence to: lestup@bgc-jena.mpg.de, linamaesu@gmail.com
# #### GitHub repository: https://github.com/linamaes/Regional_ESDL
# 
# This script does the following:
# - Load biotic unit map from NetCDF
# - Use cubefittable for calculating mean and standard deviation of the seasonal cycle
# - Compute the seasonality ratio between the annual and semiannual oscillation from the Fast Fourier power spectrum 
# - Compute the fraction of annual and semiannual oscillation from the Fast Fourier power spectrum
# 
# About the notebook:
# - It is written in Julia 1.3
# - "#" comments in the code are intended to explain specific aspects of the coding
# - New steps in workflows are introduced with bold headers
# 
# April 2021, Max Planck Insitute for Biogeochemistry, Jena, Germany

using ESDL

using ESDLPlots

using CSV

using Dates

using NPZ

using NetCDF

using Dates

using Plots

using Distributed

using OnlineStats

using WeightedOnlineStats

using FFTW

gr(size=(600,400))
default(fmt = :png)

cin = loadCube("pcacomstdTS1pc1km2014_mask")

# Note: MSC dates are loss when saving the Cube
cmsc = loadCube("pcacom1msc1km2014_mask")

# Note: SDSC dates are loss when saving the Cube
csdsc = loadCube("pcacom1sdsc1km2014_mask")

# Note: MSC dates are loss when saving the Cube
cmscfft = loadCube("pcafftmsc1kmstd2014_mask")

# Note: SDSC dates are loss when saving the Cube
csdscfft = loadCube("pcafftsdsc1kmstd2014_mask")

plotMAP(cmscfft, MSC=-29219.0)

pathin = "/my_path/.../"

path = string(pathin,"ESDLreg/netdcf/UB_IAvH_wgs84.nc");

layer = ncread(path, "layer");

cbu = CubeMem(CubeAxis[getAxis("Lon", cmsc),getAxis("Lat", cmsc)], layer)

plotMAP(cbu, dmin=0, dmax=70)

lat = (4,6)
lon = (-73,-71)
t1 = (2001:2005)

csub = cin[lat = lat, lon = lon]
busub = cbu[lat = lat, lon = lon]

plotMAP(busub)

cmsc

cbu

cmscsub = getMSC(csub)

MSCdates = cmscsub.axes[1]

indims = InDims("MSC")
outdims = OutDims(MSCdates)

function mscTime(xout,xin)
    xout[:] = xin
end

cmscyr = mapCube(mscTime, cmsc, indims=indims, outdims=outdims)

plotMAP(cmscyr, MSC=Date(1900,1,1))

tab = CubeTable(veg=cmscyr, biome=cbu, include_axes=("lat","MSC"))
meanbybiome = cubefittable(tab,WeightedMean,:veg,by=(:biome, i->(i.MSC)), weight=i->cosd(i.lat))

tab = CubeTable(veg=cmscfft, biome=cbu, include_axes=("lat","MSC"))
meanbybiomefft = cubefittable(tab,WeightedMean,:veg,by=(:biome, i->(i.MSC)), weight=i->cosd(i.lat))

tax = collect(cmscyr.axes[1].values);

meanbybiome.properties["labels"]=include("bioticunits_name.jl");

plot(tax, meanbybiome[3,:])

plot(tax, meanbybiomefft[3,:])

tab = CubeTable(veg=cmscyr, biome=cbu, include_axes=("lat","MSC"))
varbybiome = cubefittable(tab, WeightedVariance, :veg, by=(:biome, i->(i.MSC)), weight=i->cosd(i.lat))

sdbybiome = map(x -> sqrt(x), varbybiome)

tab = CubeTable(veg=cmscfft, biome=cbu, include_axes=("lat","MSC"))
varbybiomefft = cubefittable(tab, WeightedVariance, :veg, by=(:biome, i->(i.MSC)), weight=i->cosd(i.lat))

sdbybiomefft = map(x -> sqrt(x), varbybiomefft)

sdbybiome.properties["labels"]=include("bioticunits_name.jl");

### This function computes the ratio between the 2nd & (3rd+4yh) components
function spratioFx(xin)
    any(ismissing, xin) && return missing
    xfft = fft(convert(Vector{Float64}, xin))
    xout = abs(xfft[2])/(abs(xfft[3])+abs(xfft[4]))
end

ratiosp = mapslices(spratioFx, meanbybiome, dims="Category2")

# Rows are the biotic units sort by ratiosp's cube axis, and columnts are the MSC dates
meanbybiome[:,:]

function perAnnSemFx(xin)
    any(ismissing, xin) && return missing
    xfft = fft(convert(Vector{Float64}, xin))
    absx = map(x->abs(x), xfft)
    #xout = (abs(xfft[2])+abs(xfft[3])+abs(xfft[4]))/sum(absx[2:23])
    xout = sum(absx[2:4])/sum(absx[2:23])
    return xout
end

percannsem = mapslices(perAnnSemFx, meanbybiome, dims="Category2")

plot(meanbybiome[1,:])

# Concatenate data
df1 = hcat(meanbybiome.axes[1].values, ratiosp[:], percannsem[:], meanbybiome[:,:]);

# saveCube(meanbybiome, "biomepcamsc1km2014_mask")

# saveCube(meanbybiomefft, "biomepcamscfft1km2014_mask")

# saveCube(sdbybiome, "biomepcasdsc1km2014_mask")

# saveCube(sdbybiomefft, "biomepcasdscfft1km2014_mask")

# saveCube(ratiosp, "ratio_biome_pca1km2014_mask")

# saveCube(percannsem, "fraction_biome_pca1km2014_mask")
