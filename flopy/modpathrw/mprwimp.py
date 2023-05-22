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
        The model object (of type :class: ModpathRW) to which
        this package will be added.
    defaultboundary : int
        The default behavior of particles arriving to cell faces
        without connections. Allowed values are: 
         * 0: boundaries are considered open and particles are removed
         * 1: boundaries are considered closed and particles rebound
    impformat : int
        The format in which impermeable cells are defined. Allowed values are: 
         * 0: Follow flow-model inactive cells 
         * 1: Follow time-dependent flow-model inactive cells (include dry cells)
         * 2: Read icbound parameter and write with u3d 
    icbound : list, np.ndarray [int]
        Specification for impermeable cells. Particles rebound 
        in those cells with value 1. Read with u3d. 
    """
   
    
    @staticmethod
    def _ftype():
        return 'IMP'


    def __init__(
            self,
            model,
            impformat       = 0,
            defaultboundary = 0,
            icbound         = 0, # for u3d
            extension       = 'imp',
        ):
  
        # Call super constructor
        ftype      = self._ftype()
        unitnumber = model.next_unit()
        super().__init__(model, extension, ftype, unitnumber)
      
        # Default boundary 
        if ( defaultboundary not in [0,1] ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid defaultboundary {str(defaultboundary)}"
                f". Allowed values are 0 (open boundaries) "
                f"or 1 (impermeable boundaries)."
            )
        self.defaultboundary = defaultboundary

        # Impermeable cells format 
        if ( impformat not in [0,1,2] ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid impformat {str(impformat)}."
                f" Allowed values are 0 (follow ibound), "
                f"1 (follow iboundts) or 2 (given in icbound param)."
            )
        self.impformat = impformat

        # The case where icbound is needed 
        if self.impformat == 2:

            if icbound is None:
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Impermeable cells format was " 
                    f"specified as 2, but no icbound was given."
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

            if not isinstance( icbound, (int,list,np.ndarray) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for icbound. It should be list, or np.ndarray. "
                    f"{str(type(icbound))} was given."
                )
            if isinstance( icbound, int ): 
                icbound = [ icbound ]
            icbound  = np.array(icbound)
            uicbound = np.unique(icbound)
            for uc in uicbound: 
                if uc not in [0,1]:
                    raise ValueError( 
                        f"{self.__class__.__name__}:"
                        f" Invalid values for icbound specification. "
                        f"It should only contain values 0 or 1. "
                        f"The value {str(uc)} was found."
                    )

            try:
                self.icbound= Util3d(
                    model,
                    shape3d,
                    np.int32,
                    icbound,
                    name="ICBOUND",
                    locat=self.unit_number[0],
                )
            except:
                raise Exception(
                    f"{self.__class__.__name__}:"
                    f" Error while initializing distributed variable icbound. "
                    f"Is the input shape consistent with flow model dimensions ?"
                )

        # Add package 
        self.parent.add_package(self)

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
        with open(self.fn_path, 'w') as f:

            # Write the default boundary 
            f.write(f"{self.defaultboundary}\n")  

            # Write the impermeable cells format 
            f.write(f"{self.impformat}\n")

            # For impformat 2, write icbound
            if ( self.impformat == 2):
                f.write(self.icbound.get_file_entry())

        # Done
        return
