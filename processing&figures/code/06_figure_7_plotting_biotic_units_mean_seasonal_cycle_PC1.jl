# ## Main figures
# 
# ## Figure 7: Plots of the Mean Seasonal Cycle (MSC) of the first principal component by biotic units (bu)
# 
# ### Estupinan-Suarez, et al. (2021). A Regional Earth System Data Lab for Understanding Ecosystem Dynamics: An Example from Tropical South America. Front. Earth Sci. 9:613395. doi: 10.3389/feart.2021.613395
# #### Correspondence to: lestup@bgc-jena.mpg.de, linamaesu@gmail.com
# #### GitHub repository: https://github.com/linamaes/Regional_ESDL
# 
# 
# Steps:
# - Load computed MSC (mean and std) for 1st PCA component
# - Load ratio annual and semiannual ocsillation cube by bu
# - Load percentage of variance explained by bu
# - Plots MSC by biotic units. It corresponds to Figure 7A to 7D 
# 
# 
# About the notebook:
# - It is written in Julia 1.3
# - "#" comments in the code are intended to explain specific aspects of the coding
# - New steps in workflows are introduced with bold headers
# 
# April 2021, Max Planck Insitute for Biogeochemistry, Jena, Germany

using OnlineStats

using ESDL

using CSV

using Dates

using NPZ

using NetCDF

using Plots

using ESDLPlots

using DelimitedFiles

gr(size=(600,400))
default(fmt = :png)

pathin = "/my_path_in/.../"

percannsem = loadCube("fraction_biome_pca1km2014_mask")

ratiosp = loadCube("ratio_biome_pca1km2014_mask")

meanbybiome = loadCube("biomepcamsc1km2014_mask")

sdbybiome = loadCube("biomepcasdsc1km2014_mask")

cin = loadCube("pcacomstdTS1pc1km2014_mask")

lcbybuin = readdlm(string(pathin, "dataout/lc2014bybupercentage.csv"), ',');

# Convert missings from text to data missings
lcbybu = map(x->x=="missing" ? missing : x, lcbybuin)

# Exclude missings and NaNs
lcbybu2 = lcbybu[3:end,:];

cin = getMSC(cin[lat=(4,6), lon=(-73,-71), time=2001:2005])

taxin = collect(cin.axes[1].values);

tax = map(x->Dates.dayofyear(x), taxin[1:46]);

bunames = collect(ratiosp.axes[1].values);

meanbybiome

ratiosp

# Dictionary without Spanish tilde
dict = include("bioticunits_name_notilde.jl");
dictlc2 = include("lc_legend_abr2.jl");

# Map dictionary to biotic units
buid = collect(meanbybiome.axes[1].values)
buname = map(x->dict[x], buid);

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

# Smooth MSC from biotic unit for plotting purpose
smsc = mapslices(smoothMSC,meanbybiome,dims="Category2",ns=3)

# Geta axis Biotic Units ID
mscbuax = collect(meanbybiome.axes[1].values);

row = 2

buloc = 6
buid = Integer(mscbuax[buloc])

pathoutts = string(pathin,"/plots/bumscpcastd/mask/2014/axesequalsmoothLC/")
varname = "Seasonal oscillation (PCA comp. No. 1)"
nameout = "PCAvegstd"
buidx = 5
ylim1 = -3
ylim2 = 4

titlen = dict[buid]

px = plot(tax, (smsc[:,buloc]*-1),
    color="green", legend=false,
    ribbon=sdbybiome[buloc,:], fillalpha=0.1,
    lw=3, xlabel="Day of year", xtickfont=font(7),
    ylabel=varname, ytickfont=font(7), guidefont=font(8),
    title=titlen, titlefontsize=11, legendfontsize=7,
    size=(300,300), ylim=(ylim1,ylim2), dpi=200)

nameout2 = join(map(x -> isspace(titlen[x]) ? "" : titlen[x], 1:length(titlen)))

# Initialize px5 and px28
px5, px28 = 0, 0

for buloc = 2:size(ratiosp.axes[1])[1]

    buid = Integer(mscbuax[buloc])

#   Select the two most dominant land cover classes
    # Find  the row loction in lcbybu for selected by biotic unit "buidx"
    lcall = findall(row->row==buid, lcbybu2[:,1])
    # Filter land covers for the selected biotic unit
    lcsub = hcat(lcbybu2[lcall,2:end]' |> Matrix, lcbybu[1,2:end])
    # Sort land cover classes ascendingly including missings
    lcmis = sortslices(lcsub, dims=1, rev=true)
    # Exclude land cover classes with missing values
    idmis = findmax(findall(x->ismissing(x), lcmis[:,1]))[1]
    lcmax = lcmis[idmis+1:end,:]

      titlex = string(dictlc2[lcmax[1,2]], " (", round(lcmax[1,1], digits=2), "%)
  & ", dictlc2[lcmax[2,2]], " (", round(lcmax[2,1], digits=2), "%)")

      coordout = string(titlex)

      titlen = dict[buid]

#   Smoothing function
     smsc = mapslices(smoothMSC, meanbybiome, dims="Category2",ns=3)

     px = plot(tax, meanbybiome[buloc,:].*-1, color="orange",
         lw=1.5, fillalpha=0.1, legend=:false)

#   Add values of the ratio and fraction to the plot
    annotate!([
            #(180,-4.3,text(titlen,7,:center,:gray)), # -> activate to add biotic units name on the bottom
           (300,3.2,text(string("\nRatio  = " , round(ratiosp[buloc], digits=2),
            "\nFraction = ", round(percannsem[buloc], digits=2)),9))])

    plot!(tax, (smsc[:,buloc]*-1),
    color="green", legend=false,
    ribbon=sdbybiome[buloc,:], fillalpha=0.1,
    lw=3, xlabel="Day of year", xtickfont=font(7),
    ylabel=varname, ytickfont=font(7), guidefont=font(8),
    title=coordout, titlefontsize = 11, legendfontsize = 7,
    size=(300,300), ylim=(ylim1,ylim2), dpi=200)

#   Delete blank spaces from biotic units names
    nameout2 = join(map(x -> isspace(titlen[x]) ? "" : titlen[x], 1:length(titlen)))

#     savefig(px, string(pathoutts, nameout, "_bu_", buid, "_", nameout2, ".png"))


#   Save two examples for demostrative purpose
#     buid == 28 ? px28=px : buid == 5 ? px5=px : print(0)
    if buid == 5
        px5=px
        println("example ", buid)
    elseif buid == 28
        px28=px
        println("example ", buid)
    end


end

titlen = dict[buid]

buid, nameout2

px

px5

px28
