## flopyrw
An extension of [FloPy](https://github.com/modflowpy/flopy) to write input simulation files for [MODPATH-RW](https://gitub.com/upc-ghs/modpath-rw) with Python.

### Overview
Provides classes extended from the `modpath` module in `flopy` adapted to specific structures required by MODPATH-RW. Also introduces new package writers required by the program, consistent with the [Documentation of Input-Output](https://github.com/upc-ghs/modpath-rw/doc/modpath-rw_IO_v100_.pdf). 


### Quickstart
**Install**

To install the package from source, clone the repository:

```
git clone https://github.com/upc-ghs/flopyrw
```
and install 

```
pip install -e /the/path/to/flopyrw/
```

You can also install the current release from [PyPI](https://pypi.org/project/flopyrw/):

```
pip install flopyrw
```

**Use it**

Classes follow the same logic than `flopy`, configuring packages on top of a MODFLOW flow-model object.

```py
from flopyrw import modpathrw

# Suppose an existing gwt (mf6) model with 
# the aux variable 'CONCENTRATION' in 'WEL-1' package... 

# Create a modpathrw model
mprw = modpathrw.ModpathRW(modelname="mprwsim", flowmodel=gwt) 

# Random walk options
modpathrw.ModpathRWOpts(
    mprw,
    dimensionsmask=[1,1,0], # Random walk in x,y and not in z
    timestep='min',
    ctdisp=0.1,
    courant=0.1
)

# Dispersion parameters 
modpathrw.ModpathRWDsp( 
    mprw,
    alphal= 0.1,
    alphat= 0.01, 
    dmeff = 0.0, 
)

# bas 
modpathrw.ModpathRWBas(mp,porosity=0.3)

# Define a solute source 
modpathrw.ModpathRWSrc(
    mp,
    sources=(
        "WEL-1", # package name
        [   # auxvar, particlesmass, release template
            [ "CONCENTRATION", 300.0, (4,4,1)], 
        ],
    ),
)

# Configure the simulation 
simconfig = {
    'simulationtype'    : 'rwendpoint', 
    'trackingdirection' : 'forward',
    'weaksinkoption'    : 'stop_at',
    'weaksourceoption'  : 'pass_through',
    'stoptimeoption'    : 'specified',
    'stoptime'          : 1.0,
}
modpathrw.ModpathRWSim(
    mprw, 
    **simconfig
)

# Write the input files
mprw.write_input()

# And run 
mprw.run_model()
```

**Note**: In order to run a model via the interface a [MODPATH-RW](https://gitub.com/upc-ghs/modpath-rw) executable is required. 


### Testing
A suite of [automated tests](autotest/) is available verifying different aspects of the interface and the program. In order to run these tests, the current release candidate of [FloPy](https://github.com/modflowpy/flopy) is required. Install it with the command:

```
pip install https://github.com/modflowpy/flopy/zipball/develop
```

Install the additional test dependencies with:

```
pip install ".[test]"
```

You can follow the [FloPy test guidelines](https://github.com/modflowpy/flopy/blob/develop/DEVELOPER.md#running-tests) for running and debugging tests. 

Run the complete test suite from the folder ``autotest`` with the command:

```
pytest -s -v 
```

## Resources
* [MODPATH](https://www.usgs.gov/software/modpath-particle-tracking-model-modflow)
* [modpath-v7 repository](https://github.com/MODFLOW-USGS/modpath-v7)
* [modpath-omp repository](https://github.com/MARSoluT/modpath-omp)
* [flopy](https://github.com/modflowpy/flopy)
