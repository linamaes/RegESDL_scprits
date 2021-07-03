#!/usr/bin/env python
# coding: utf-8

# ## Data access
#
# ## Script 02d01: Code to download the regional ESDL using Python
#
#
# ### Estupinan-Suarez, Gans, et al. (2021). A Regional Earth System Data Lab for Understanding Ecosystem Dynamics: An Example from Tropical South America. Front. Earth Sci. 9:613395. doi: 10.3389/feart.2021.613395
# #### Correspondence to: lestup@bgc-jena.mpg.de, linamaesu@gmail.com
# #### GitHub repository: https://github.com/linamaes/Regional_ESDL
#
#
# About the script:
# - This code was written in jupyter notebooks
# - This notebook exemplifies how to access the RegESDL. Here, we showcase how to download the RegESDL using s3fs in Python for local use
#
# About the notebook:
# - The notebook kernel is Python 3.7.1, and it uses the s3fs package version 0.5.2
# - "#" comments in the code are intended to explain specific aspects of the coding
# - New steps in workflows are introduced with bold headers
# - For complete information of the datasets check the orignal data source documentation. Citations are available on the referred manuscript and supplementary tables.
#
# March 2021, Max Planck Insitute for Biogeochemistry, Jena, Germany

# ## Import requiered packages


import s3fs


# ## Commands used to download the RegESDL

# In[ ]:


# Create a handle to the otc object store
s3 = s3fs.S3FileSystem(
    anon = True,
)


# In[ ]:


# List of available colombia data cubes:
s3.ls('esdl-esdc-v2.0.1')


# In[ ]:


# Commands used to download a cube
obs.get(obs.ls('esdl-esdc-v2.0.1/Cube_2019highColombiaCube_184x120x120.zarr'),'./mylocalhighrescube', recursive=True);


# Commands to access the data cube using xcube and then continue analysis using xarray:

# In[ ]:


import xarray as xr
import xcube


# In[ ]:


from xcube.core.dsio import open_dataset

ds = open_dataset("https://s3.eu-central-1.amazonaws.com/esdl-esdc-v2.0.1/Cube_2019highColombiaCube_184x120x120.zarr", s3_kwargs=dict(anon=True))


# In[ ]:
