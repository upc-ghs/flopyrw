## flopyrwpt
`python` interface for writing input files and postprocessing tools for transport models based on MODPATH-RW


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


## Writing MODPATH-RW input files

Package extends from the `flopy` classes for writing `modpath` input files. For `modpath-rw`, additional classes are considered for writing dispersion, reconstruction and other configuration files required for the program. 

## Resources

* [MODPATH](https://www.usgs.gov/software/modpath-particle-tracking-model-modflow)
* [modpath-v7 repository](https://github.com/MODFLOW-USGS/modpath-v7)
* [modpath-omp repository](https://github.com/MARSoluT/modpath-omp)
* [flopy](https://github.com/modflowpy/flopy)
