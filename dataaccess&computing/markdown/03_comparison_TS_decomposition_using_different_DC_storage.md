## Computational Performance

## Script 03: Computational comparison of time series decomposition using different Data Cubes storage

### Estupinan-Suarez, Gans, et al. (2021). A Regional Earth System Data Lab for Understanding Ecosystem Dynamics: An Example from Tropical South America. Front. Earth Sci. 9:613395. doi: 10.3389/feart.2021.613395
#### Correspondence to: lestup@bgc-jena.mpg.de, linamaesu@gmail.com
#### GitHub repository: https://github.com/linamaes/Regional_ESDL


This script illustrates the time difference on time series computing when using different cube storage versions. We provide two examples:
    
- Example 1: Time series are decomposed in different time-scales i.e. Longer-term, seasonal cycle and short(fast) oscillation, plus the trend 
- Example 2: It finds the maximum value along a time period at pixel level



About the notebook
- It is written in Julia 1.3
- "#" comments in the code are intended to explain specific aspects of the coding
- New steps in the workflow are introduced with bold headers

March 2021, Max Planck Insitute for Biogeochemistry, Jena, Germany


## Load packages


```julia
using ESDL
```


```julia
using ESDLPlots
```


<script>
// Immediately-invoked-function-expression to avoid global variables.
(function() {
    var warning_div = document.getElementById("webio-warning-1067746827823222194");
    var hide = function () {
        var script = document.getElementById("webio-setup-16373082021116835020");
        var parent = script && script.parentElement;
        var grandparent = parent && parent.parentElement;
        if (grandparent) {
            grandparent.style.display = "none";
        }
        warning_div.style.display = "none";
    };
    if (typeof Jupyter !== "undefined") {
        console.log("WebIO detected Jupyter notebook environment.");
        // Jupyter notebook.
        var extensions = (
            Jupyter
            && Jupyter.notebook.config.data
            && Jupyter.notebook.config.data.load_extensions
        );
        if (extensions && extensions["webio-jupyter-notebook"]) {
            // Extension already loaded.
            console.log("Jupyter WebIO nbextension detected; not loading ad-hoc.");
            hide();
            return;
        }
    } else if (window.location.pathname.includes("/lab")) {
        // Guessing JupyterLa
        console.log("Jupyter Lab detected; make sure the @webio/jupyter-lab-provider labextension is installed.");
        hide();
        return;
    }
})();

</script>
<p
    id="webio-warning-1067746827823222194"
    class="output_text output_stderr"
    style="padding: 1em; font-weight: bold;"
>
    Unable to load WebIO. Please make sure WebIO works for your Jupyter client.
    For troubleshooting, please see <a href="https://juliagizmos.github.io/WebIO.jl/latest/providers/ijulia/">
    the WebIO/IJulia documentation</a>.
    <!-- TODO: link to installation docs. -->
</p>



## Set reading and writing paths


```julia
pathin = "../mypath/"# => Set your path
```




    "../mypath/"



## Load Cubes


```julia
# Cube storage version for temporal analysis
ctem = Cube(string(pathin,"/Cube_2020highColombiaCube_184x120x120.zarr/"))
```




    Collection of ZArray Cube with the following dimensions
    Lon                 Axis with 2760 Elements from -82.99622135 to -60.00464665
    Lat                 Axis with 3360 Elements from 13.99613735 to -13.99541735
    Time                Axis with 782 Elements from 2001-01-05T00:00:00 to 2017-12-31T00:00:00
    Variable            Axis with 5 elements: lai fpar gross_primary_productivity evi ndvi 
    Total size: 168.85 GB





```julia
# Cube storage version for spatial analysis
cspa = Cube(string(pathin,"/Cube_2020highColombiaCube_1x3360x2760.zarr/"))
```




    Collection of ZArray Cube with the following dimensions
    Lon                 Axis with 2760 Elements from -82.99622135 to -60.00464665
    Lat                 Axis with 3360 Elements from 13.99613735 to -13.99541735
    Time                Axis with 782 Elements from 2001-01-05T00:00:00 to 2017-12-31T00:00:00
    Variable            Axis with 5 elements: lai fpar gross_primary_productivity evi ndvi 
    Total size: 168.85 GB




## Settings for cube subset


```julia
yrstart = 2001
yrend = 2005
latsub =(0,1)
```




    (0, 1)



## Example 1: Time series decomposition


```julia
@time ctemfft = filterTSFFT(ctem[time=yrstart:yrend, lat=latsub])
```

    [32mProgress: 100%|█████████████████████████████████████████| Time: 0:02:55[39m


    186.832809 seconds (1.57 G allocations: 55.861 GiB, 6.66% gc time)





    Collection of ZArray Cube with the following dimensions
    Time                Axis with 230 Elements from 2001-01-05T00:00:00 to 2005-12-31T00:00:00
    Scale               Axis with 4 elements: Trend Long-Term Variability Annual Cycle Fast Oscillations 
    Lon                 Axis with 2760 Elements from -82.99622135 to -60.00464665
    Lat                 Axis with 120 Elements from 0.9961893499999996 to 0.004526649999999694
    Variable            Axis with 5 elements: lai fpar gross_primary_productivity evi ndvi 
    Total size: 7.09 GB





```julia
@time cspafft = filterTSFFT(cspa[time=yrstart:yrend, lat=latsub])
```

    ┌ Warning: There are still cache misses
    └ @ ESDL.DAT /Net/Groups/BGI/scratch/lestup/julia_atacama_depots/packages/ESDL/skMpG/src/DAT/DAT.jl:608
    ┌ Warning: There are compressed caches misses, you may want to use a different cube chunking
    └ @ ESDL.DAT /Net/Groups/BGI/scratch/lestup/julia_atacama_depots/packages/ESDL/skMpG/src/DAT/DAT.jl:609
    [32mProgress: 100%|█████████████████████████████████████████| Time: 0:36:33[39m


    2194.006997 seconds (1.53 G allocations: 929.078 GiB, 4.11% gc time)





    Collection of ZArray Cube with the following dimensions
    Time                Axis with 230 Elements from 2001-01-05T00:00:00 to 2005-12-31T00:00:00
    Scale               Axis with 4 elements: Trend Long-Term Variability Annual Cycle Fast Oscillations 
    Lon                 Axis with 2760 Elements from -82.99622135 to -60.00464665
    Lat                 Axis with 120 Elements from 0.9961893499999996 to 0.004526649999999694
    Variable            Axis with 5 elements: lai fpar gross_primary_productivity evi ndvi 
    Total size: 7.09 GB




## Example 2: Maximum values along time series


```julia
@time ctempax= mapslices(maximum, ctem[time=yrstart:yrend, lat=latsub], dims="time")
```

    [32mProgress: 100%|█████████████████████████████████████████| Time: 0:00:41[39m


     43.233203 seconds (756.82 M allocations: 17.970 GiB, 9.68% gc time)





    In-Memory data cube with the following dimensions
    Lon                 Axis with 2760 Elements from -82.99622135 to -60.00464665
    Lat                 Axis with 120 Elements from 0.9961893499999996 to 0.004526649999999694
    Variable            Axis with 5 elements: lai fpar gross_primary_productivity evi ndvi 
    Total size: 7.9 MB





```julia
@time cspamax = mapslices(maximum, cspa[time=yrstart:yrend, lat=latsub], dims="time")
```

    ┌ Warning: There are still cache misses
    └ @ ESDL.DAT /Net/Groups/BGI/scratch/lestup/julia_atacama_depots/packages/ESDL/skMpG/src/DAT/DAT.jl:608
    ┌ Warning: There are compressed caches misses, you may want to use a different cube chunking
    └ @ ESDL.DAT /Net/Groups/BGI/scratch/lestup/julia_atacama_depots/packages/ESDL/skMpG/src/DAT/DAT.jl:609
    [32mProgress: 100%|█████████████████████████████████████████| Time: 0:12:38[39m


    758.984782 seconds (754.09 M allocations: 249.939 GiB, 4.20% gc time)





    In-Memory data cube with the following dimensions
    Lon                 Axis with 2760 Elements from -82.99622135 to -60.00464665
    Lat                 Axis with 120 Elements from 0.9961893499999996 to 0.004526649999999694
    Variable            Axis with 5 elements: lai fpar gross_primary_productivity evi ndvi 
    Total size: 7.9 MB





```julia

```
