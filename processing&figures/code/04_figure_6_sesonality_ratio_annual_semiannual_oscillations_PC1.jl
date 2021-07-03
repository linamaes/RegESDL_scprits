# ## Main figures
# 
# ## Figure 6: (6B) Seasonality ratio map of annual and semiannual oscillations from the Mean Seasonal Cycle (MSC) and (6A) plots from three pixels
# 
# ### Estupinan-Suarez, et al. (2021). A Regional Earth System Data Lab for Understanding Ecosystem Dynamics: An Example from Tropical South America. Front. Earth Sci. 9:613395. doi: 10.3389/feart.2021.613395
# #### Correspondence to: lestup@bgc-jena.mpg.de, linamaesu@gmail.com
# #### GitHub repository: https://github.com/linamaes/Regional_ESDL
# 
# This script does the following:
# - Computes the Mean Seasonal Cycle (MSC) of the first principal component from the PCA analysis. For more details about PCA analysis see script 03 - PCA processing.
# - Imports and applies a mask based on MODIS quality flags. See the mask generation script 01.
# - Calcultes the ratio between annual and semiannual oscillation based on the Fast Fourier Spectrum.
# 
# About the notebook
# - It is written in Julia 1.3
# - "#" comments in the code are intended to explain specific aspects of the coding
# - New steps in workflows are introduced with bold headers
# 
# March 2021, Max Planck Institute for Biogeochemistry, Jena, Germany

using ESDL

using ESDLPlots

using Plots

using FFTW

using Distributed

using CSV

using DelimitedFiles

gr(size=(600,400))
default(fmt = :png)

pathmask = ".../my_pathin/"
pathout = ".../my_pathout/"

pcacom = loadCube("pcacomstdTS1pc1km2014_mask") # Component 1

plotMAP(pcacom, time=Date(2007,12,31), dmin=-5, dmax=5)

cmsc = getMSC(pcacom)

plotMAP(cmsc, MSC=Date(1900,1,1), dmin=-5, dmax=5)

# Load mask
maskin = CSV.read(string(pathmask, "msc_missings_mask.csv"), header=false);

# Convert NaNs to missings and transpose matrix
mask = map(x-> isnan(x) ? missing : x==0 ? missing : x, Matrix(maskin));
mask = mask' |> Array;

function maskout(xout, xin, mask)
    xout[:] = xin.*mask
end

indims = InDims("Lon","Lat")
outdims = OutDims("Lon","Lat")

cmscm = mapCube(maskout, cmsc, mask, indims=indims, outdims=outdims)

plotMAP(cmscm, MSC=Date(1900,1,1), dmin=-5, dmax=5)

# This fuction computes the ratio between 2nd and (3rd+4th) fft frequency power
# from fftw (fft[2]/(fft[3]+fft[4])

function spratioFx2(xin)
    any(ismissing, xin) && return missing
    xfft = fft(convert(Vector{Float64}, xin))
    xout = abs(xfft[2])/(abs(xfft[3])+abs(xfft[4]))
end

ratiosp2 = mapslices(spratioFx2, cmscm, dims="MSC")

# Colors differ with the publication colorscheme. The publication map was done using matplotlib in Python
plotMAP(ratiosp2, dmin=0, dmax=2)

ratioout2 = map(x -> ismissing(x) ? NaN : x, (ratiosp2[:,:]' |> Matrix));
# writedlm(string(pathout,"ratio_annual_semiannual_maskosc234.csv"), ratioout2, ",")

function perAnnSemFx(xin)
    any(ismissing, xin) && return missing
    xfft = fft(convert(Vector{Float64}, xin))
    absx = map(x->abs(x), xfft)
    xout = (abs(xfft[2])+abs(xfft[3])+abs(xfft[4]))/sum(absx[2:23])
    return xout
end

percentannsem = mapslices(perAnnSemFx, cmsc, dims="MSC")

plotMAP(percentannsem)

addprocs(8);

@everywhere using ESDL

cfft = filterTSFFT(pcacom)

cfftmsc = getMSC(cfft[Scale="Annual Cycle"])

rmprocs(workers())

cfftmscm = mapCube(maskout, cfftmsc, mask, indims=indims, outdims=outdims)

plotMAP(cfftmscm, MSC=Date(1900,1,1))

addprocs(6);

pcacom

@everywhere function getsdSC(xout, xin, timeax)
    any(ismissing, xin) && return missing
    timedoy = map(x->Dates.dayofyear(x), timeax)
    #xout = zeros(46)
    for i in 1:46
        idx = findall(x-> x==timedoy[i], timedoy)
        xout[i] = std(map(x->xin[x], idx))
    end
    return xout
end

axmsc = cmsc.axes[1]

indimssdsc = InDims("Time")
outdimssdsc = OutDims(axmsc)

@everywhere using Dates, OnlineStats

tall = collect(pcacom.axes[1].values);

csdsc = mapCube(getsdSC, pcacom, tall, indims=indimssdsc, outdims=outdimssdsc)

csdscfft = mapCube(getsdSC, cfft[Scale="Annual Cycle"], tall, indims=indimssdsc, outdims=outdimssdsc)

 rmprocs(workers())

csdscm = mapCube(maskout, csdsc, mask, indims=indims, outdims=outdims)

csdscfftm = mapCube(maskout, csdscfft, mask, indims=indims, outdims=outdims)

workers()

plotMAP(csdscm, MSC=Date(1900,1,1))

plotMAP(csdscfftm, MSC=Date(1900,12, 27))

# saveCube(ratiosp2, "ratio_annual_semiannual_pca1kmstd2014_mask")

# saveCube(cmscm, "pcacom1msc1km2014_mask")

# saveCube(cfftmscm, "pcafftmsc1kmstd2014_mask")

# saveCube(csdscm, "pcacom1sdsc1km2014_mask")

# saveCube(csdscfftm, "pcafftsdsc1kmstd2014_mask")

# Get the longitude and latitude coordinates
lonaxc = collect(cmsc.axes[2].values)
lataxc = collect(cmsc.axes[3].values);

# Define days of the year from MSC
taxin = collect(cmsc.axes[1].values);
tax = map(x->Dates.dayofyear(x), taxin);

# foldp = "plots/pxs"

# Define function for plotting

function plotpx(clon, clat)

    titlen = "\nRatio = "#"\nRatio (ann/sem) = "
    titlen2 = "  Fraction = " #"\nPercentage (ann&sem) = "


    coordout = string("Px. coord: Lat. ", round(lataxc[clat], digits=2) ,"°, Lon. ", round(lonaxc[clon], digits=2),"°",
        titlen, round(ratiosp2[clon,clat], digits=2),
        titlen2, round(percentannsem[clon,clat], digits=2))

    p1 = plot(tax, cmsc[:,clon, clat].*-1, label  = "PC1 MSC", color="black",
        ribbon=csdsc[:,clon,clat], fillalpha=0.1, lw=1,
        xlabel="Day of year", xguidefontsize=14, xtickfontsize=15,
        ylim=(-3, 3), yguidefontsize=14,  ytickfontsize=15, #ylabel="PC1",
        title=coordout, titlefontsize = 12, legendfontsize = 7, legend=:false,
        size=(300, 300)) #ylim=(3, 8.52

    plot!(tax, cfftmsc[:, clon, clat].*-1, label = "FFT (annual)", lw=3, color="orange")

#     savefig(p1, string(rootp,foldp,"lon_",clon,"lat_",clat,"2014v2_mask",dpix,"dpi.png"))

    return p1
end

px1 = plotpx(1500, 1000)

px2 = plotpx(1570, 1550)

px3 = plotpx(1560, 1700)

px4 = plotpx(1300,1750)
