'''
Configuration of MODPATH-RW basic package
'''

# python 
import numpy as np

# flopy 
from flopy.modpath import Modpath7Bas

class ModpathRWBas( Modpath7Bas ):
    '''
    MODPATH-RW Basic class 

    Extends flopy.modpath.Modpath7Bas
    '''


    @staticmethod
    def _ftype():
        return 'MPRWBAS'


    def __init__(self, *args, **kwargs):
        from copy import copy

        # Add a layer of validation for porosity values
        if 'porosity' in kwargs.keys():
            if ( not isinstance( kwargs['porosity'], (int,float,list,np.ndarray)) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for porosity. It should be real values, but "
                    f"{str(type(kwargs['porosity']))} was given."
                )
            # take porosity just for this validation
            # and pass the kwarg to the super constructor
            porosity = copy(kwargs['porosity'])
            if ( isinstance( porosity, (int,float) ) ):
                porosity = [porosity]
            if ( isinstance( porosity, list ) ): 
                porosity = np.array(porosity)
            if ( porosity.dtype not in [np.float32, np.float64, np.int32, np.int64 ] ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for porosity. It should be real values, but "
                    f"{str(porosity.dtype)} was given."
                )
            # Verify that all values are within a valid range
            if ( np.any( porosity < 0.0 ) or np.any( porosity > 1.0 ) ): 
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid values for porosity. It should contain real values " 
                    f"only in the range [0.0,1.0]."
                )
            # Complain if only zeros
            if ( np.all( porosity == 0.0 ) ): 
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid values for porosity. It should contain real values " 
                    f"only in the range [0.0,1.0], but only zeros were found."
                )

        # Call parent constructor
        super().__init__(*args,**kwargs)

        # Not the most elegant solution but consistent
        ftype = self._ftype()
        self.name  = [ftype]
        self._name = [ftype]
        self._generate_heading()
