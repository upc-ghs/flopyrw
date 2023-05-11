'''
Configuration of MODPATH-RW species
'''

# python
import numpy as np

# flopy
from flopy.pakbase import Package

# local 
from .utils import count_instances # Increment COUNTER
from .mprwdsp import ModpathRWDsp


class ModpathRWSpc( Package ):
    """
    MODPATH-RW Species Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    dispersion: ModpathRWDsp object
        Defines dispersion model and properties to this specie.
    pgroups: list, np.array
        Particle groups related to this specie.
    """


    # Class properties
    COUNTER    = 0
    UNITNUMBER = 0
    INSTANCES  = []


    @count_instances
    def __init__(
        self,
        model,
        dispersion, # it's like a foreign key
        pgroups    = None,
        id         = None,
        stringid   = None,
        extension  = 'spc',
    ):
        
        # Define UNITNUMBER if the first instance is created
        if self.__class__.COUNTER == 1:
            unitnumber = model.next_unit()
            super().__init__(model, extension, "SPC", unitnumber)
            self.__class__.UNITNUMBER = self.unit_number[0]
        else:
            # Needed for Util3d
            self._parent = self.INSTANCES[0]._parent

        # Assign dispersion
        if not isinstance( dispersion, ModpathRWDsp ):
            raise Exception('flopyrw:ModpathRWSpc: dispersion parameter should be of ModpathRWDsp type. Given ', type( dispersion ) )
        self.dispersion = dispersion

        # Assign pgroups id's
        if pgroups is None: 
            warnings.warn("flopyrw:ModpathRWSpc: no pgroups were specified, defaults to pgroups = [0]")
            pgroups = [0]
            # Or leave as None
        if not isinstance(pgroups, (list,np.array)):
            raise Exception('flopyrw:ModpathRWSpc: pgroups should be of type list or np.array. Given ', type( pgroups ) )
        if isinstance( pgroups, list ):
            pgroups = np.array( pgroups )
        # Requires some more sanity checks to verify integer id's and so on
        # Might be worth considering removing particle groups altogether and extend 
        # species from those classes.
        self.pgroups = pgroups
        

        # Define id
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER

        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = 'SPC'+str(self.__class__.COUNTER)


        # Add package and save instance        
        if self.__class__.COUNTER == 1:
            self.parent.add_package(self)
        self.__class__.INSTANCES.append( self )


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

        # Open file for writing
        f = open(self.INSTANCES[0].fn_path, "w")

        # Write how many INSTANCES 
        # and loop over them
        f.write(f"{self.COUNTER}\n")

        for ins in self.INSTANCES:

            # Write id
            f.write(f"{ins.id}\n")

            # Write stringid
            f.write(f"{ins.stringid}\n")

            # Need to write the related groups 
            # only if particlesmassoption not equal 2
            # If particlesmassoption equal 2, then 
            # soluteId is infered from given particlegroups
            if ( self.INSTANCES[0]._parent.particlesmassoption != 2 ):
                # Write the number of groups
                f.write(f"{len(ins.pgroups)}\n")

                for idpg, pg in enumerate(ins.pgroups):
                    # Write the pgroup id in the list of pgroups
                    f.write(f"{pg+1}\n")

            # Write dispersion "foreign key"        
            f.write(f"{ins.dispersion.id}\n")
            f.write(f"{ins.dispersion.stringid}\n")


        # Done
        return
