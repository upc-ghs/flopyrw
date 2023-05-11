'''
Configuration of MODPATH-RW impermeable cells
'''

# python
import numpy as np

# flopy
from flopy.pakbase import Package
from flopy.utils import Util3d


class ModpathRWImp( Package ):
    """
    MODPATH-RW Impermeable Cells Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    defaultboundary : int
        The default behavior of particles arriving to cell faces without connections.
        Allowed values are: 
            0: boundaries are considered open and particles arriving to inactive cells are removed
            1: boundaries are considered closed and particles rebound
    impformat : int
        The format in which impermeable cells are defined. Allowed values are: 
            0: Follow flow-model inactive cells 
            1: Follow time-dependent flow-model inactive cells (include dry cells)
            2: Read icbound parameter and write with u3d 
    icbound : int, np.ndarray
        Specification for impermeable cells. Particles rebound in those cells with value 1. 
    """
    
    def __init__(
            self,
            model,
            impformat       = 0,
            defaultboundary = 0,
            icbound         = 0, # for u3d
            extension       = 'imp',
        ):
    
        # Define UNITNUMBER if the first instance is created
        unitnumber = model.next_unit()
        super().__init__(model, extension, "IMP", unitnumber)
      
        # Default boundary 
        if ( defaultboundary not in [0,1] ): 
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid defaultboundary ' +str(defaultboundary)+
                '. Allowed values are 0 (open boundaries) or 1 (impermeable boundaries).'
            )
        self.defaultboundary = defaultboundary

        # Impermeable cells format 
        if ( impformat not in [0,1,2] ): 
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid impformat ' +str(impformat)+
                '. Allowed values are 0 (follow ibound), 1 (follow iboundts) or 2 (given in icbound param).'
            )
        self.impformat = impformat

        if self.impformat == 2: 
            if icbound is None:
                raise ValueError(
                    self.__class__.__name__ + ':' + 
                    ' Impermeable cells format was specified as 2, but no icbound was given.' 
                )

            shape = model.shape
            if len(shape) == 3:
                shape3d = shape
            elif len(shape) == 2:
                shape3d = (shape[0], 1, shape[1])
            else:
                shape3d = (1, 1, shape[0])
            self.model = model
            self.shape3d = shape3d

            self.icbound= Util3d(
                model,
                shape3d,
                np.int32,
                icbound,
                name="ICBOUND",
                locat=self.unit_number[0],
            )

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
        f = open(self.fn_path, "w")

        # Write the default boundary 
        f.write(f"{self.defaultboundary}\n")  

        # Write the impermeable cells format 
        f.write(f"{self.impformat}\n")

        # For impformat 2, write icbound
        if ( self.impformat == 2):
            f.write(self.icbound.get_file_entry())

        # And close
        f.close()

        return
