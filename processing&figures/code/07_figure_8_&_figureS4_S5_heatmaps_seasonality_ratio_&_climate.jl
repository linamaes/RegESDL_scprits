# ## Main and supplementary figures
# 
# ## Figure 7, S4, S5: Heatmaps of seasonality ratio of annual/semiannual oscillations and climate variables
# 
# 
# ### Estupinan-Suarez, et al. (2021). A Regional Earth System Data Lab for Understanding Ecosystem Dynamics: An Example from Tropical South America. Front. Earth Sci. 9:613395. doi: 10.3389/feart.2021.613395
# #### Correspondence to: lestup@bgc-jena.mpg.de, linamaesu@gmail.com
# #### GitHub repository: https://github.com/linamaes/Regional_ESDL
# 
# This script does the following:
# - Loads cubes with seasonality ratio (annual/semiannual oscillation), biotic units map, land cover data by biotic units, and climate variables (precipitation of the driest month, maximum temperature of the warmest month, median annual cloud frecuency)
# - Select climate data by biotic units
# - Plot heatmaps of climate variables versus seasonality ratio
# 
# About the notebook:
# - It is written in Julia 1.3
# - "#" comments in the code are intended to explain specific aspects of the coding
# - New steps in workflows are introduced with bold headers
# 
# 
# April 2021, Max Planck Institute for Biogeochemistry, Jena, Germany

using OnlineStats

using ESDL

using ESDLPlots

using Plots

using NetCDF

using DelimitedFiles

gr(size=(400,350)) # gr(size=(600,400)) #
default(fmt = :png)

# Location RegESDL
pathin1 = "/my_pathin1/.../"

#Location land cover data by biotic units
pathin2 = "/my_pathin2/.../"

# Location Biotic units map in NetCDF format
pathin3 = "/my_pathin3/.../"

lat = (4,6)
lon = (-73,-71)

time = (2001:2015)

c = Cube(string(pathin1,"/Cube_2019highColombiaCube_184x120x120.zarr/"))

# foreach(println,caxes(c)[4].values)

# Percentage of land cover by biotic untis
lcbybuin = readdlm(string(pathin2,"dimred/dataout/lc2014bybupercentage.csv"), ',')
lcbybu = map(x->x=="missing" ? missing : x, lcbybuin)
lcbybu2 = lcbybu[3:end,:];

dictlc2 = include("lc_legend_abr2.jl");

cmask = c[variable="fapar", Time=Date(2001,1,5)]

mask1 = map(x->ismissing(x) ? missing : 1, cmask)

layer = ncread(string(pathin3,"/UB_IAvH_wgs84.nc"), "layer");

cbu = CubeMem(CubeAxis[getAxis("Lon", c),getAxis("Lat", c)], layer)
cbu.properties["labels"]=include("bioticunits_name_notilde.jl");

plotMAP(cbu, dmin=0, dmax=70)

bunames = include("bioticunits_name_notilde.jl")
i = 1
bunames[i]

# This is a descriptive variable. It only has data on the first time-step
cppt = c[variable="PrecipitationofDriestMonth", Time=Date(2001,1,5)]

plotMAP(cppt)

ctem = c[variable="MaxTemperatureofWarmestMonth", Time=Date(2001,1,5)]

heatmap(ctem[:,:]'[end:-1:1,:], aspect_ratio = :equal)

ccloud = c[variable=["MODCF_monthlymean_01","MODCF_monthlymean_02","MODCF_monthlymean_03",
        "MODCF_monthlymean_04","MODCF_monthlymean_05","MODCF_monthlymean_06",
        "MODCF_monthlymean_07","MODCF_monthlymean_08","MODCF_monthlymean_09",
        "MODCF_monthlymean_11","MODCF_monthlymean_10","MODCF_monthlymean_12"],
    Time=Date(2001,1,5)]#, Time=Date(2001,1,5)]

# Calculate annual median of clouds coverage
medianax = CategoricalAxis("Median", ["percentage"])
indims = InDims("Variable")
outdims = OutDims(medianax)

function medianFx(xout, xin)
    #@show size(xin)#, size(xin)
    sort!(xin)
    xout .= (xin[6]+xin[7])/200 # coefficient 2  * 100 to get %
end

ccloudmd = mapCube(medianFx, ccloud, indims=indims, outdims=outdims)

# Plot data
heatmap(ccloudmd[1,:,:]'[end:-1:1,:], aspect_ratio = :equal)

mask1

ratin = readdlm(string(pathin2, "/monobimodal/dataout/ratio_annual_semiannual_maskosc234.csv"),',')
ratmis = map(x->isnan(x) ? missing : x, ratin' |> Array)
rat1 = CubeMem(CubeAxis[getAxis("Lon", c),getAxis("Lat", c)], ratmis)

# rat1 = loadCube("ratio_annual_semiannual_pca1kmstd2014")

plotMAP(rat1, dmin=0, dmax=2)

# Define a mask for vegetation variables, ratio map, and climate variable
function maskvarratioFx(mask1, rat1, cin)
    maskall = map((x,y,z)->x*y*z, mask1, rat1, cin)
    mask2 = map(x->ismissing(x) ? missing : 1, maskall[:,:]);
end

# Apply mask to ratio map based on the assessed variable
mask2ppt = maskvarratioFx(mask1, rat1, cppt)
mask2tem = maskvarratioFx(mask1, rat1, ctem)
mask2cloud = maskvarratioFx(mask1, rat1, ccloudmd[Median="percentage"]);

heatmap(mask2ppt[:,:]'[end:-1:1,:], aspect_ratio = :equal)

indims = InDims("Lon","Lat")
outdims = OutDims("Lon","Lat")

function maskout(xout, xin, mask)
    xout[:] = xin.*mask
end

# Exclude pixels from the ratio map that do not match the assessed variable using the respective mask
# ratm = mapCube(maskout, rat1, mask2, indims=indims, outdims=outdims)

ratmppt = mapCube(maskout, rat1, mask2ppt, indims=indims, outdims=outdims)
ratmtem = mapCube(maskout, rat1, mask2tem, indims=indims, outdims=outdims)
ratmclouds = mapCube(maskout, rat1, mask2cloud, indims=indims, outdims=outdims)

heatmap(ratmppt[:,:]'[end:-1:1,:], aspect_ratio = :equal)

# Exclude pixels from the assessed variable that do not match the ratio map
# varm = mapCube(maskout, cvar, mask2, indims=indims, outdims=outdims)

varmppt = mapCube(maskout, cppt, mask2ppt, indims=indims, outdims=outdims)
varmtem = mapCube(maskout, ctem, mask2tem, indims=indims, outdims=outdims)
varmcloud = mapCube(maskout, ccloudmd, mask2cloud, indims=indims, outdims=outdims)

heatmap(varmppt[:,:]'[end:-1:1,:], aspect_ratio = :equal)

cbu

i = 27
m1 = map(x -> x==i ? 1.0 : 0.0, cbu)

plotMAP(m1)

idx = findall(x-> x==i, cbu[:,:])
cbuin = cbu[:,:]
validx = (map(x->cbuin[x], idx));

rat1d = ratmppt[:,:];

varin = varmppt[:,:];
varname = "Precipitation driest month (mm)"

# Number of pixels for selected biotic unit
size(collect(skipmissing(map(x->rat1d[x], idx))))

ratvaridx = hcat((collect(skipmissing(map(x->rat1d[x], idx)))), (collect(skipmissing(map(x->varin[x], idx)))))
#collect(skipmissing(ratidx))
ratvaridx[:,1] = map(x-> x>2 ? 2 : x, ratvaridx[:,1]);

titleout = size(ratvaridx)[1]

histogram2d(ratvaridx[:,1], ratvaridx[:,2], title = string(bunames[i], " \n(n = ", size(ratvaridx)[1] , ")"),
    yaxis=(varname), xaxis=("Seasonality ratio (annual/semiannual)"), c=:blues)

function plothmbu(bunames, cbu, rat1d, varin, varname, pathout, yst, yend)

    # From 1 until the length of biotic units names
    for i = 1:66

        titlen = bunames[i]

        # Find the two most dominant land cover clases for the selected biotic unit
        lcall = findall(row->row==i, lcbybu2[:,1])
        lcsub = hcat(lcbybu2[lcall,2:end]' |> Matrix, lcbybu[1,2:end])
        lcmis = sortslices(lcsub, dims=1, rev=true)
        idmis = findmax(findall(x->ismissing(x), lcmis[:,1]))[1]
        lcmax = lcmis[idmis+1:end,:]

        # Set figure title based on dominant land covers
        titlex = string(dictlc2[lcmax[1,2]], " (", round(lcmax[1,1], digits=2), "%) +
",dictlc2[lcmax[2,2]], " (", round(lcmax[2,1], digits=2), "%)")

        # Select pixels of interest bases on the selected biotic unit
        idx = findall(x-> x==i, cbu[:,:])
        ratvaridx = hcat((collect(skipmissing(map(x->rat1d[x], idx)))),
            (collect(skipmissing(map(x->varin[x], idx)))))

        if minimum(ratvaridx[:,2]) < 0
            @show i
        end

    # Plot commands
        pout = histogram2d(ratvaridx[:,1], ratvaridx[:,2],
            titlefontsize = 13, title = string(titlex, "  (n=", size(ratvaridx)[1] , ")"),
            xlim =(0, 3), xaxis=("Seasonality ratio (annual/semiannual)"),  xtickfont=font(11),
            yaxis=(varname), ytickfont=font(12, "TimesNewRoman"),
#              ylim = (yst, yend), #activate for equal axis range e.g. temperature and clouds
            xguidefont=font(15, "TimesNewRoman"), yguidefont=15,
            dpi=200)

        # Save figure
        nameout2 = join(map(x -> isspace(titlen[x]) ? "" : titlen[x], 1:length(titlen)))
        savefig(pout, string(pathout, "heatmap_bu_", i, "_", nameout2, ".png"))

    end

end

# Variable name for the y-axis label
varname = "Precip. driest month (mm)"

# Location for saving figures
locfol = "minpremonth/axisdiffer/"
pathout = string(pathin2, "dimred/plots/heatmaps/mask/2014/", locfol);

varname

plothmbu(bunames, cbu, ratmppt[:,:], varmppt[:,:], varname, pathout, 0, 800)

# Variable name for the y-axis label
varname = "Max. °T warmest month (°C)"

# Location for saving figures
locfol = "maxtemwarmestmon/axisequal/original/"
pathout = string(pathin2, "dimred/plots/heatmaps/mask/2014/", locfol);

plothmbu(bunames, cbu, ratmtem[:,:], varmtem[:,:], varname, pathout, 0, 35)

# Variable name for the y-axis label
varname = "Median cloud frecuency (%)"

# Location for saving figures
locfol = "clouddays/median/axisequal/original/"
pathout = string(pathin2, "dimred/plots/heatmaps/mask/2014/", locfol);

plothmbu(bunames, cbu, ratmclouds[:,:], varmcloud[:,:,1], varname, pathout, 40, 100)
