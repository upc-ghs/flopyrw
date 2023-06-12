## flopyrw
Interface for writing input simulation files for MODPATH-RW, based on [FloPy](https://github.com/modflowpy/flopy)

### Overview
Provides classes extended from the `modpath` module at `flopy` and further functionalities to write the package input files for [MODPATH-RW](https://gitub.com/upc-ghs/modpath-rw).

Classes follow the same logic than `flopy`, configuring packages based on a MODFLOW flow-model object.

### Quicktstart

Install the package 

Use it

```py
from flopy import modpathrw

# Create a modpathrw model

# Append some packages

# Configure the simulation 

# Write the files

# And run 
```

Note: the interface requires the MODPATH-RW executable. 






## Install

Clone reposistory and install as dependency in your environment. Suggested process is 

```
git clone https://gitlab.com/upc-ghs/flopyrw.git
```

Change directory to the package

```
cd flopyrw
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
pip install -e /the/path/to/flopyrw/
```

Then, published changes in the central repository can be integrated by regularly checking the output of the command (executed inside the cloned `flopyrw` repository)

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
