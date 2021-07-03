# ## Supplementary figures
# 
# ## Figure S3: Loadings of principal components
# 
# ### Estupinan-Suarez, et al. (2021). A Regional Earth System Data Lab for Understanding Ecosystem Dynamics: An Example from Tropical South America. Front. Earth Sci. 9:613395. doi: 10.3389/feart.2021.613395
# #### Correspondence to: lestup@bgc-jena.mpg.de, linamaesu@gmail.com
# #### GitHub repository: https://github.com/linamaes/Regional_ESDL
# 
# This scripts does the following:
# - Plots loadings from the first principal component. Colors and map format differe from the puclication, because the paper figure was produced using matplotlib in Python
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

pcasiglod = loadCube("pcaComSigLoadings1km2014_mask")

abs(-8)

# Get absoulte values for pc1
pc1 = map(x->abs(x), pcasiglod[info_PCA="LoadingPC1"])

# GPP
plotMAP(pc1, Components_or_Var="sigPC_var1", colorm=colormap("greens"), dmin=0, dmax=0.8)

# NDVI
plotMAP(pc1, Components_or_Var="sigPC_var2", colorm=colormap("greens"), dmin=0, dmax=0.8)

# EVI
plotMAP(pc1, Components_or_Var="sigPC_var3", colorm=colormap("greens"), dmin=0, dmax=0.8)

# FPAR
plotMAP(pc1, Components_or_Var="sigPC_var4", colorm=colormap("greens"), dmin=0, dmax=0.8)
