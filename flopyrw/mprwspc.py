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
        The model object (of type :class: ModpathRW) to which
        this package will be added.
    dispersion: ModpathRWDsp object
        Defines dispersion model and properties for this solute. Is a required foreign key. 
    pgroups : list, np.array (int)
        Particle groups related to this specie. These are interpreted only if 
        particlesmassoption != 2. These ids are zero based (while writing is corrected to one based)
        and are related to the order in which the particle groups are specified to the program.
        In case the simulation is configured with particlesmassoption == 2, 
        then the particle groups associated to a specie are interpreted (internally)
        from the specified solute ids for each particle group.
    stringid : str
        String identifier for this species. 
    extension : str
        The extension for the species file, by default spc. 
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
        id         = None, # internal
        stringid   = None,
        extension  = 'spc',
    ):
        
        # Define UNITNUMBER if the first instance is created
        if self.__class__.COUNTER == 1:
            unitnumber = model.next_unit()
            super().__init__(model, extension, "SPC", unitnumber)
            self.__class__.UNITNUMBER = self.unit_number[0]
        else:
            self._parent = self.INSTANCES[0]._parent

        # Assign dispersion
        if not isinstance( dispersion, ModpathRWDsp ):
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid dispersion parameter of type ' +str(type(dispersion))+
                '. It should be type ModpathRWDsp. '
            )
        self.dispersion = dispersion

        # Will be written only if particlesmassoption != 2
        if pgroups is not None:
            if not isinstance(pgroups, (list,np.array)):
                raise ValueError(
                    self.__class__.__name__ + ':' + 
                    ' Invalid pgroups parameter of type ' +str(type(pgroups))+
                    '. It should be of type list or np.array.'
                )
            if isinstance( pgroups, list ):
                pgroups = np.array( pgroups )
        # Requires some more sanity checks to verify integer id's.
        # Might be worth considering a verification versus already 
        # specified particle groups, which in the end would require
        # forcing giving this package after IC, SRC.
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

            # Write stringid
            f.write(f"{ins.stringid}\n")

            # Need to write the related groups 
            # only if particlesmassoption not equal 2
            # If particlesmassoption equal 2, then 
            # soluteId is infered from given particlegroups
            if ( self.INSTANCES[0]._parent.particlesmassoption != 2 ):
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
                        self.__class__.__name__ + ':' + 
                        ' Invalid pgroups parameter of type ' +str(type(ins.pgroups))+
                        '. It should be of type list or np.array with at least one element.'
                    )

            # Write dispersion "foreign key", only if
            # speciesdispersionoption == 1
            if  self.INSTANCES[0]._parent.speciesdispersionoption == 1:
                f.write(f"{ins.dispersion.stringid}\n")


        # Done
        return
