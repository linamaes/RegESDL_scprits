# ## Main figures
#
# ## Figure 4: Variance explained by principal components of vegetation variables
#
# ### Estupinan-Suarez, et al. (2021). A Regional Earth System Data Lab for Understanding Ecosystem Dynamics: An Example from Tropical South America. Front. Earth Sci. 9:613395. doi: 10.3389/feart.2021.613395
# #### Correspondence to: lestup@bgc-jena.mpg.de, linamaesu@gmail.com
# #### GitHub repository: https://github.com/linamaes/Regional_ESDL
#
# This scipts does the following:
# - Gap fills and standarized data to zero mean and variance of 1
# - Computes Principal Component Analysis (PCA) for vegetation variables i.e., GPP, NDVI, EVI, FAPAR
# - Calculates the significance of each principal component (PC)
# - Extracts loading values of each component
# - Computes the mean seasonal cycle (MSC) and standard deviation seasonal cycle (SDSC) of the first PC
#
# About the notebook:
# - It is written in Julia 1.3
# - "#" comments in the code are intended to explain specific aspects of the coding
# - New steps in workflows are introduced with bold headers
#
# April 2021, Max Planck Institute for Biogeochemistry, Jena, Germany

using ESDL

using ESDLPlots

using CSV

using Missings

using Plots

using MultivariateStats

using OnlineStats

using DelimitedFiles

using Distributed

using NPZ

using Statistics

gr(size=(600,400))
default(fmt = :png)

pathin = "/my_path/.../"

# High spatial resolution
c = Cube(string(pathin, "Cube_2020highColombiaCube_184x120x120.zarr/"))

# Low spatial resolution
# c = Cube(string(pathin, "Cube_2020lowColombiaCube_184x60x60.zarr/"))

csub = c[time=2001:2014, var=["gross_primary_productivity", "ndvi",
        "evi","fpar"]]#

# Set time for plotting
t = Date(2014, 1, 31)

plotMAP(csub, time=t, variable="gross_primary_productivity", dmin=0, dmax=10)

pathcsv = "/my_mask_path/.../"

pathcsv = "/Net/Groups/BGI/people/lestup/ESDLreg/output/csv/"

maskin = CSV.read(string(pathcsv, "msc_missings_mask.csv"), header=false);

maskin = Matrix(maskin);

mask = map(x-> isnan(x) ? missing : x==0 ? missing : x, maskin' |> Array);

csub

# lonaxis = RangeAxis("Lon", lon)
# lataxis = RangeAxis("Lat", lat)
lonaxis = csub.cubeaxes[1]
lataxis = csub.cubeaxes[2]

# Create a cube from the mask
cmask = CubeMem(CubeAxis[lonaxis, lataxis], mask)

plotMAP(cmask)

workers()

## add workers
addprocs(6)

csub

csub2 = map(x-> ismissing(x) ? missing : x==255.0 ? missing : x==65335.0 ? missing : x, csub)

#csub2 = mapCube(maskout, csub[lat=lat, lon=lon], mask, indims=indims, outdims=outdims)

plotMAP(csub2, time=t, variable="gross_primary_productivity", dmin=0, dmax=10)

#cfill = gapFillMSC(csub)
@time cfill = gapFillMSC(csub2)

cfill

plotMAP(cfill, time=t, variable="gross_primary_productivity", dmin=0, dmax=10)

@everywhere using OnlineStats

# Variables are standardized to zero mean and variance of 1

@everywhere function stdFx(xout, xin)
    xout[:] .= map(x -> (x - OnlineStats.mean(xin))/sqrt(var(xin)), xin)
end

indimsstd = InDims("Time")
outdimsstd = OutDims("Time")

#@time
cstd = mapCube(stdFx, cfill, indims = indimsstd, outdims = outdimsstd)

plotMAP(cstd, time=t, variable="gross_primary_productivity")

# addprocs(8)

workers();

# rmprocs(workers())

@everywhere using MultivariateStats

@everywhere using ESDL

@everywhere function pcaSigFx(xout, xin0::Any, loopvars, maxoutdim=4, method=:svd)
    #@show loopvars
    if any(ismissing.(xin0))
        return missing
    else
       #@show size(xin0)
       xin = convert(Array{Float32}, xin0' |> Matrix)
        # projection(pca) give the rotation matrix
        xin_pca = fit(PCA, xin; maxoutdim=maxoutdim,  method=method, pratio=1.0)
       #@show size(xout)
        xout[:,1] .= (principalvars(xin_pca) ./ tprincipalvar(xin_pca))
        xout[:,2:4] .= projection(xin_pca)[:,1:3]
        return xout

    end
end

sigaxis = CategoricalAxis("info_PCA", ["PC_Significance","LoadingPC1", "LoadingPC2",
        "LoadingPC3"])#,, "LoadingPC4"

varsigaxis = CategoricalAxis("Components_or_Var", ["sigPC_var1","sigPC_var2","sigPC_var3", "sigPC_var4"]) #,"sigPC_var4"

indims = InDims("Time","Variable")
outdims = OutDims(varsigaxis, sigaxis)

cstd

########################################
# NOTE for axes interpretation:
# 'Components_or_Var' has a different meaning based on the 'info_PCA' class, as follows:
#  For 'PC_Significance' the values are  significance of PC1, PC2, PC3 and PC4 respectively
#  For 'LoadingsPC1/2/3' the values are the loadings of variable 1, 2, 3, and 4 respectively
#  Check the variables order from the input data.
########################################

@time pcasiglod = mapCube(pcaSigFx, cstd, 4, indims=indims, outdims=outdims)

# saveCube(pcasiglod, "pcaComSigLoadings1km2014_cube2020")

indims = InDims("Lon","Lat")
outdims = OutDims("Lon","Lat")

rmprocs(workers())

indims = InDims("Lon","Lat")
outdims = OutDims("Lon","Lat")

# This function assigns values of missing to low quality pixels

function maskout(xout, xin, mask)
    xout[:] = xin.*mask
end

pcasiglodout = mapCube(maskout, pcasiglod, mask, indims=indims, outdims=outdims)

# saveCube(pcasiglodout, "pcaComSigLoadings1km2014_mask")

pcasigm = pcasiglodout[info_PCA="PC_Significance"]

plotMAP(pcasigm[Components="sigPC_var1"], colorm=colormap("reds"), dmin=0, dmax=1,
    oceancol=colorant"white",
    misscol=colorant"white",
    dpi=150)

plotMAP(pcasigm[Components="sigPC_var2"], colorm=colormap("greens"), dmin=0, dmax=1,
    oceancol=colorant"white",
    misscol=colorant"white")

plotMAP(pcasigm[Components="sigPC_var3"], colorm=colormap("blues"), dmin=0, dmax=1,
    oceancol=colorant"white",
    misscol=colorant"white") # misscol=colorant"white",

minimum(skipmissing(pcasigm[:,:,1])), maximum(skipmissing(pcasigm[:,:,1]))

minimum(skipmissing(pcasigm[:,:,2])), maximum(skipmissing(pcasigm[:,:,2]))

minimum(skipmissing(pcasigm[:,:,3])), maximum(skipmissing(pcasigm[:,:,3]))

quantile(skipmissing(pcasigm[:,:,1]),[0.01, .99])

quantile(skipmissing(pcasigm[:,:,2]),[0.01, .99])

quantile(skipmissing(pcasigm[:,:,3]),[0.01, .91])

pcasigm = pcasiglodout[info_PCA="PC_significance"]

p1 = plotMAPRGB(pcasigm, rgbaxis="Components",
    oceancol=colorant"white",
    misscol=colorant"white",
    c1="sigPC_var1",
    c2="sigPC_var2",
    c3="sigPC_var3",
    cType=RGB)

addprocs(6);

workers();

@everywhere using MultivariateStats

@everywhere function pcaFx(xout, xin0::Any, maxoutdim=maxd, method=:svd)
    if any(ismissing.(xin0))
        return missing
    else
       xin = convert(Array{Float32}, xin0' |> Matrix)
        # projection(pca) give the rotation matrix
        xin_pca = fit(PCA, xin; maxoutdim=maxoutdim,  method=method, pratio=1.0)
        t = transform(xin_pca, xin)
        xout .= missing
        for i in 1:min(size(xout,2),size(t,1))
            xout[:,i] .= t[i,:]
        end
    end
end


axispca = CategoricalAxis("Components", ["C1","C2","C3","C4"]) #,"C5"

indimspca = InDims("Time","Variable")
outdimspca = OutDims("Time",axispca)

cstd

#pcacom = mapCube(pcaFx, cstd, 5, indims=indims, outdims=outdims)
@time pcacom = mapCube(pcaFx, cstd, 4, indims=indimspca, outdims=outdimspca)

plotMAP(pcacom, time=t, Component="C1")

# saveCube(pcacom[Components="C1"],"pcacomstdTS1pc1km2014_mask")

@time cmsc = getMSC(pcacom[Components="C1"])

@everywhere function maskout(xout, xin, mask)
    xout[:] = xin.*mask
end

indims = InDims("Lon","Lat")
outdims = OutDims("Lon","Lat")

workers();

rmprocs(workers())

@time cmscm = mapCube(maskout, cmsc, mask, indims=indims, outdims=outdims)

plotMAP(cmscm, MSC=Date(1900, 12, 27))

# saveCube(cmscm, "pcacom1msc1km2014_mask")

axmsc = cmsc.axes[1]

@everywhere function getSDSC(xout, xin, timeax)
    any(ismissing, xin) && return missing
    timedoy = map(x->Dates.dayofyear(x), timeax)
    #xout = zeros(46)
    for i in 1:46
        idx = findall(x-> x==timedoy[i], timedoy)
        xout[i] = std(map(x->xin[x], idx))
    end
    return xout
end

indimssdsc = InDims("Time")
outdimssdsc = OutDims(axmsc)

@everywhere using Dates, OnlineStats

tall = collect(pcacom.axes[1].values);

pcacom[Components="C1"]

@time csdsc = mapCube(getSDSC, pcacom[Components="C1"], tall, indims=indimssdsc, outdims=outdimssdsc)

rmprocs(workers())

@time csdscm = mapCube(maskout, csdsc, mask, indims=indims, outdims=outdims)

plotMAP(csdscm, MSC=Date(1900, 12, 27))

# saveCube(csdscm,"pcacom1sdsc1km2014_mask")
