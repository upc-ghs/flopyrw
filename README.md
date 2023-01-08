# flopyrwpt
`python` interface for writing input files for solute transport models based on MODPATH-RW, based in `flopy`.

## Overview
Package provides classes similar to those found in the `flopy/modpath` module for writing input files. These include
    
* `ModpathRwptDispersion`
* `ModpathRwptReconstruction`
* `ModpathRwptObs`
* `ModpathRwptSim`

Similar to `flopy`, these classes group parameters and input file writing logic for the different model options.

## Install

Clone reposistory and install as dependency in your environment. Suggested process is 

```
git clone https://gitlab.com/upc-ghs/flopyrwpt.git
```

Change directory to the package

```
cd flopyrwpt
```

Activate the `python` environment

```
source /path/to/env/bin/activate
```

And install

```
python setup.py install
```

Another alternative would be to directly install the cloned repository once the `python` environment is activated. This is done with the command

```
pip install -e /the/path/to/flopyrwpt/
```

Then, published changes in the central repository can be integrated by regularly checking the output of the command (executed inside the cloned `flopyrwpt` repository)

```
git fetch 
```

In case of changes, these can be integrated with `git pull`.


## Writing MODPATH-RW input files
Package extends from the `flopy` classes for writing `modpath` input files. For `modpath-rw`, additional classes are considered for writing dispersion, reconstruction and other configuration files required for the program. Workflow for applying these packages follows the same logic than a `modpath-v7` simulation in `flopy`. 



## Resources

* [MODPATH](https://www.usgs.gov/software/modpath-particle-tracking-model-modflow)
* [modpath-v7 repository](https://github.com/MODFLOW-USGS/modpath-v7)
* [modpath-omp repository](https://github.com/MARSoluT/modpath-omp)
* [flopy](https://github.com/modflowpy/flopy)
