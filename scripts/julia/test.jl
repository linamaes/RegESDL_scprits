using ESDL

using ESDLPlots

using Dates

using OnlineStats

using DelimitedFiles

using Distributed

pathin = "../mypath/"# => Set your path

pathout = "../mypath"# => Set your path


# High resolution
c = Cube(string(pathin,"/Cube_2020highColombiaCube_184x120x120.zarr/"))


# Set assessed period based on variables temporal overlapping
c = c[time=2001:2014]

# Set time for plotting
t = Date(2014, 1, 31)

# Visualization of data subset
plotMAP(c, time=t, variable="gross_primary_productivity", dmin=0, dmax=10)



mask9999 = map(x->ismissing(x) ? 9999 : x, maskall[:,:]')

# writedlm(string(pathout,"msc_missings_mask.csv"), maskout, ",")

# writedlm(string(pathout,"msc_missings_mask9999.csv"), mask9999, ",")

