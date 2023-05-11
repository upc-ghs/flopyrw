'''
Configuration of MODPATH-RW initial conditions
'''

# python
import numpy as np

# flopy
from flopy.pakbase import Package
from flopy.utils import Util3d

# local 
from .utils import count_instances # Increment COUNTER


class ModpathRWIc( Package ):
    """
    MODPATH-RW Initial Condition Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    concentration : ndarray
        Interpreted as resident concentration, that is, the product beween porosity and
        dissolved concentration. 
    """


    # Class properties
    COUNTER    = 0
    UNITNUMBER = 0
    INSTANCES  = []


    @count_instances
    def __init__(
        self,
        model,
        kind          = 1,
        mass          = .1,
        concentration = None,
        soluteid      = 1, # replace by specie ?
        id            = None,
        stringid      = None,
        extension     = 'ic',
    ):
       

        # Define UNITNUMBER if the first instance is created
        if self.__class__.COUNTER == 1:
            unitnumber = model.next_unit()
            super().__init__(model, extension, "IC", unitnumber)
            self.__class__.UNITNUMBER = self.unit_number[0]
        else:
            # Needed for Util3d
            self._parent = self.INSTANCES[0]._parent


        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])
        self.model = model
        self.shape3d = shape3d

        # Determines how to write the initial condition
        if (kind not in [1]): 
            raise ValueError('flopyrw:ModpathRWIc: Invalid value for kind, should be 1. Given ', str(kind) )
        self.kind = kind

        # Need health checks
        self.mass       = mass
        self.soluteid   = soluteid

        # define id
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER

        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = 'IC'+str(self.__class__.COUNTER)

        # and concentration: needs id
        if (
            ( concentration is not None ) and 
            ( self.kind == 1 )
        ):
            if not isinstance( concentration, np.ndarray ):
                raise TypeError('flopyrw:ModpathRWIc: concentration should be an np.ndarray')
            self.concentration = Util3d(
                model,
                shape3d,
                np.float32,
                concentration,
                name=self.stringid,
                locat=self.__class__.UNITNUMBER,
            )
        elif ( (self.kind ==1) and (concentration is None) ):
            raise ValueError('flopyrw:ModpathRWIc: kind 1 requires concentration array, None given.')
     

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

            # Write string id 
            f.write(f"{ins.stringid}\n")

            # Kind/format of initial condition 
            f.write(f"{ins.kind}\n")

            # 1: resident concentration array  
            if ins.kind == 1:

                # Give particles mass 
                f.write(f"{ins.mass:.10f}\n")

                # If solutes are being specified, do it
                if( 
                    ( self.INSTANCES[0]._parent.particlesmassoption == 2 ) or
                    ( self.INSTANCES[0]._parent.speciesdispersionoption == 1 ) 
                ):
                    f.write(f"{ins.soluteid}\n")

                # And concentration distribution
                f.write(ins.concentration.get_file_entry())


        # Done
        return
