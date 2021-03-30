using ESDL

using ESDLPlots

using Dates

using OnlineStats

using DelimitedFiles

using Distributed

pathin = "../mypath/"# => Set your path


################### pathin = "/Net/Groups/BGI/work_2/ColombiaCube/"

# pathout = "../mypath"# => Set your path
################### 
pathout = "/Net/Groups/BGI/people/lestup/ESDLreg/output/csv/"

# High resolution
c = Cube(string(pathin,"/Cube_2020highColombiaCube_184x120x120.zarr/"))


# Set assessed period based on variables temporal overlapping
c = c[time=2001:2014]

# Set time for plotting
t = Date(2014, 1, 31)

# Visualization of data subset
plotMAP(c, time=t, variable="gross_primary_productivity", dmin=0, dmax=10)

# Because NDVI and EVI have the same quality flags, we just calculate a general mask for both
# vegetation indices (vi). The missing value from MODIS vegetation indices is "65335.0"
cvi = map(x-> ismissing(x) ? missing : x==65335.0 ? missing : x, c[variable="ndvi"])

# The missing value from fpar is "255.0"
cfpar = map(x-> ismissing(x) ? missing : x==255.0 ? missing : x, c[variable="fpar"])

# Get the number of time steps per year
npy = getNpY(cvi)

# Define cube dimmensions for applying the MapCube function 

# Set input dimension
indims = InDims("time")

# Define "Day-of-year (DOY) as the new output dimension"
newaxis = RangeAxis("DOY",collect(1:npy))

# Set output dimension
outdims = OutDims(newaxis)

# @everywhere 
# This function calculates the number of missing values for each MSC time-step 
function doymisFx(xout, xin, idoy, doyin)# idoy, doyin, 
    #size(xout)
    for i in 1:length(idoy)
        idx = findall(x->x==idoy[i], doyin)
        xout[i] = sum(ismissing.(xin[idx]))        
    end
    return xout
end

# Extract the time axes from the cube
tax = collect(cvi.axes[3].values)

# Define the day-of-year for each date along the time axis 
doyin = map(x->Dates.dayofyear(x), tax)

# Find the unique day-of-year along the time axis and sort it
idoy = sort!(unique(doyin))

# addprocs(4)
# workers()
# rmprocs(workers())

# Map function to calculate missing values for each annual time step
cvimis = mapCube(doymisFx, cvi, idoy, doyin, indims=indims, outdims=outdims)

plotMAP(cvimis, DOY=23)
# Commands for saving the output and loading it 
savecube(cviannual, "missings4msc_timesteps_vi_2001to2014")
cviannual = loadcube("missings4msc_timelsteps_vi_2001to2014")
cfparmis = mapCube(doymisFx, cfpar, idoy, doyin, indims=indims, outdims=outdims)

plotMAP(cfparmis, DOY=23)
# Commands for saving the output and loading it back
savecube(cfparannual, "missings4msc_timsteps_fpar_2001to2014")
cfparannual = loadcube("missings4msc_timesteps_fpar_2001to2014")
cgppmis = mapCube(doymisFx, c[variable="gross_primary_productivity"], 
    idoy, doyin, indims=indims, outdims=outdims)

plotMAP(cgppmis, DOY=23)

cvimismax = mapslices(maximum, cvimis, dims="DOY")

plotMAP(cvimismax)

cfparmismax = mapslices(maximum, cfparmis, dims="DOY")

plotMAP(cfparmismax)

cgppmismax = mapslices(maximum, cgppmis, dims="DOY")

plotMAP(cgppmismax)

vimask = map(x->ismissing(x) ? missing : x==14 ? 0 : 1, cvimismax)

plotMAP(vimask)

fparmask = map(x->ismissing(x) ? missing : x==14 ? 0 : 1, cfparmismax)

plotMAP(fparmask)

gppmask = map(x->ismissing(x) ? missing : x==14 ? 0 : 1, cgppmismax)

plotMAP(gppmask)

maskall = map((x,y,z)->x.*y.*x, vimask, fparmask, gppmask)

plotMAP(maskall)

maskout = map(x->ismissing(x) ? NaN : x, maskall[:,:]')

mask9999 = map(x->ismissing(x) ? 9999 : x, maskall[:,:]')

# writedlm(string(pathout,"msc_missings_mask.csv"), maskout, ",")

# writedlm(string(pathout,"msc_missings_mask9999.csv"), mask9999, ",")

pathout

workers()


