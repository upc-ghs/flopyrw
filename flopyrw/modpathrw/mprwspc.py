'''
Configuration of MODPATH-RW species
'''

# python
import numpy as np

# flopy
from flopy.pakbase import Package

# local 
from .utils import multipackage
from .mprwdsp import ModpathRWDsp


class ModpathRWSpc( Package ):
    """
    MODPATH-RW Species Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class: ModpathRW) to which
        this package will be added.
    dispersion: ModpathRWDsp object (required)
        Relates a dispersion model to this solute. Is a required foreign key. 
    pgroups : list, np.ndarray [int]
        Particle groups related to this specie. These are interpreted 
        only if particlesmassoption != 2. These ids are zero based 
        (while writing is corrected to one based) and are related to 
        the order in which the particle groups are specified to the program.
        In case the simulation is configured with particlesmassoption == 2, 
        then the particle groups associated to a specie are interpreted 
        (internally) from the specified solute ids for each particle group.
    stringid : str
        String identifier for this species. By default concatenates 
        package ftype and integer id.
    extension : str
        The extension for the species file, by default spc. 
    """


    @staticmethod
    def _ftype():
        return 'SPC'


    @multipackage
    def __init__(
        self,
        model,
        dispersion, # it's like a foreign key
        pgroups    = None,
        stringid   = None,
        extension  = 'spc',
    ):

        # Call parent constructor
        ftype = self._ftype()
        if ( model.multipackage[ftype]['count'] == 0 ): 
            super().__init__(
                model,
                extension,
                ftype,
                model.multipackage[ftype]['unitnumber'],
                allowDuplicates=True,
            )
        else:
            # Pass the parent
            self._parent = model.multipackage[ftype]['instances'][0]._parent

        # Assign dispersion
        if not isinstance( dispersion, ModpathRWDsp ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid dispersion parameter of type {str(type(dispersion))}."
                f" It should be type ModpathRWDsp."
            )
        self.dispersion = dispersion

        # Will be written only if particlesmassoption != 2
        if pgroups is not None:
            if not isinstance(pgroups, (int,list,np.ndarray)):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid pgroups parameter of type {str(type(pgroups))}."
                    f" It should be of type int, list or np.ndarray."
                )
            if isinstance( pgroups, list ):
                pgroups = np.array( pgroups )
            if isinstance( pgroups, int ):
                pgroups = np.array( [pgroups] )
            # Requires some more sanity checks to verify integer id's.
            for pg in pgroups: 
                if ( not isinstance( pg, (int,np.int32,np.int64) ) ):
                    raise TypeError(
                        f"{self.__class__.__name__}:"
                        f" Invalid pgroups parameter of type {str(type(pg))}."
                        f" It should be an integer pgroup id."
                    )

        # Might be worth considering a verification versus 
        # specified particle groups, which would require
        # checking the existence of pgroups at the simulation pkg,
        # or valid IC/SRC packages.
        self.pgroups = pgroups
      
        # String id 
        if (stringid is not None): 
            if ( not isinstance( stringid, str ) ): 
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for stringid. It should be str, but"
                    f" {str(type(stringid))} was given."
                )
            self.stringid = stringid
        else:
            self.stringid = ftype+str(self.id)

        # Done
        return


    def write_file(self, check=False):
        """
        Write the package file
        Parameters
        ----------
        check : boolean
            Check package data for common errors. (default False)
        Returns
        -------
        None
        """

        # Shortcut to multipackage dict
        mpkg = self.parent.multipackage[self._ftype()]
        
        # Write
        with open( mpkg['instances'][0].fn_path, 'w' ) as f:

            # Write how many instances...
            f.write(f"{mpkg['count']}\n")

            # ...and loop over them
            for ins in mpkg['instances']:

                # Write stringid
                f.write(f"{ins.stringid}\n")

                # Need to write the related groups 
                # only if particlesmassoption not equal 2
                # If particlesmassoption equal 2, then 
                # soluteId is infered from given particlegroups

                # BaseModel has a weird getattr method and if
                # for any reason the program is trying to write
                # a model without a defined simulation package, 
                # particlesmassoption is undefined and throws 
                # an unclear message obtained from the overloaded
                # __getattr__
                if ( not hasattr(self.parent,'particlesmassoption') ):
                    raise Exception( 
                        f"{self.__class__.__name__}:"
                        f" The particlesmassoption was not found in parent package."
                        f" Did you define a ModpathRWSim package ?"
                    )
                if ( self.parent.particlesmassoption != 2 ):
                    if ( 
                        ( ins.pgroups is not None ) and
                        ( len(ins.pgroups) > 0 )
                    ):
                        # Write the number of groups
                        f.write(f"{len(ins.pgroups)}\n")

                        for idpg, pg in enumerate(ins.pgroups):
                            # Write the pgroup id in the list of pgroups
                            # Notice the plus one.
                            f.write(f"{pg+1}\n")
                    else:
                        raise ValueError(
                            f"{self.__class__.__name__}:"
                            f" Invalid pgroups parameter of type {str(type(ins.pgroups))}."
                            f" While using particlesmassoption != 2, the ModpathRWSpc package"
                            f" requires the list of pgroups indexes related to this species."
                        )
                if ( not hasattr(self.parent,'speciesdispersionoption') ):
                    raise Exception( 
                        f"{self.__class__.__name__}:"
                        f" The speciesdispersionoption was not found in parent package." 
                        f" Did you define a ModpathRWSim package ?"
                    )
                # Write dispersion "foreign key", only if
                # speciesdispersionoption == 1
                if self.parent.speciesdispersionoption == 1:
                    f.write(f"{ins.dispersion.stringid}\n")


        # Done
        return
