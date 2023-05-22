'''
Configuration of MODPATH-RW options
'''

# python
import numpy as np

# flopy
from flopy.pakbase import Package


class ModpathRWOpts( Package ):
    """
    MODPATH-RW Options Package Class. 

    Parameters
    ----------
    model     : model object
        The model object (of type :class: ModpathRW) to which
        this package will be added.
    timestep  : str
        Define method for computing particles time step. Allowed values are:
            * 'adv'   : local time step computed with advection criteria
            * 'disp'  : local time step computed with dispersion criteria
            * 'min'   : local time step obtained as the minimum between 'adv' and 'disp'
            * 'fixed' : the same for all particles, given by the user in the 'deltat' keyword 
    courant   : float
        Courant number that will be used at each cell for computing time step with 'adv' criteria. 
        Should be greater than zero. 
    ctdisp    :  float 
        CT constant for computing time step with 'disp' criteria. Should be greater than zero.
    deltat    :  float 
        Value employed when timestep selection is 'fixed'. Should be greater than zero.
    advection : str
        Advection model to be used for advective component in RW integration.
        Could be set to 'eulerian' or 'exponential' (default is 'eulerian').
    dimensionsmask : [int,int,int]
        Determine on which dimensions RW displacements are calculated. 
        A value of 1 indicates that dimension is active and 0 inactive.
        For example, for a 2D RW model (x,y model and single layer), dimensionsmask should be [1,1,0].
        At least one dimension is required. Defaults to 3D ( [1,1,1] )
    randomgenerator : int 
        Determine the function for random number generation. Allowed values are: 
            * 0: determined by compiler (gfortran - based on random_number; ifort - Ziggurat)
            * 1: based on random_number (parallel scalability with gfortran, not with ifort)
            * 2: based on Ziggurat method (parallel scalability with gfortran and ifort)
    extension : str, optional
        File extension (default is 'rwopts').
    """


    @staticmethod
    def _ftype():
        return 'RWOPTS'


    def __init__(
        self,
        model,
        timestep         = 'min',
        courant          = 0.1,
        ctdisp           = 0.1,
        deltat           = 1.0,
        advection        = 'eulerian',
        dimensionsmask   = [1,1,1],
        randomgenerator  = 0, 
        extension        = 'rwopts',
    ):

        ftype = self._ftype()
        unitnumber = model.next_unit()
        super().__init__(model, extension, ftype, unitnumber)

        # timestep
        if ( not isinstance( timestep, str ) ): 
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for timestep selection. It should be str, but"
                f" {str(type(timestep))} was given."
            )
        if ( timestep not in ['adv', 'disp', 'min', 'fixed'] ):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Time step selection {str(timestep)} is not valid."
                f". Allowed values are: adv, disp, min, fixed."
            )
        self.timestep = timestep

        # timestep selection params
        # courant
        if ( not isinstance(courant, (int,float) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for courant parameter. It should be float/int, but"
                f" {str(type(courant))} was given."
            )
        if (courant>0): 
            self.courant  = courant
        else:
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Courant number for advection criteria should be greater than zero."
            )
        # ctdisp
        if ( not isinstance(ctdisp, (int,float) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for ctdisp parameter. It should be float/int, but"
                f" {str(type(ctdisp))} was given."
            )
        if (ctdisp>0): 
            self.ctdisp  = ctdisp
        else:
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" CT constant for dispersion criteria should be greater than zero."
            )
        # deltat
        if ( not isinstance(deltat, (int,float) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for deltat parameter. It should be float/int, but"
                f" {str(type(deltat))} was given."
            )
        if (deltat>0): 
            self.deltat  = deltat
        else:
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Time step for particles displacement should be greater than zero."
            )

        # advection 
        if ( not isinstance( advection, str ) ): 
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for advection selection. It should be str, but"
                f" {str(type(advection))} was given."
            )
        if ( advection not in ['eulerian','exponential'] ):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Advection model {str(advection)} is not valid."
                f". Allowed values are: eulerian, exponential."
            )
        self.advection  = advection

        # dimensionsmask
        if not isinstance( dimensionsmask, (list,tuple,np.ndarray) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" dimensionsmask should be of type list/tuple/np.ndarray."
                f" Given  {str(type(dimensionsmask))}"
            )
        if isinstance( dimensionsmask, (list,tuple) ):
            for dm in dimensionsmask: 
                if ( not isinstance( dm, int ) ):
                    raise TypeError(
                        f"{self.__class__.__name__}:"
                        f" Invalid type in dimensionsmask. Allowed values "
                        f"are the integers 0 or 1, but "
                        f"{str(type(dm))} was given."
                    )
            dimensionsmask = np.array(dimensionsmask).astype(np.int32)
        if ( dimensionsmask.dtype not in [np.int32,np.int64]):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type in dimensionsmask. It should be "
                f"integer, but {str(dimensionsmask.dtype)} was given."
            )
        if len(dimensionsmask) != 3:
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" dimensionsmask should be of length 3."
                f". Given {str(len(dimensionsmask))}"
            )
        for idm, dm in enumerate(dimensionsmask):
            if ( dm not in [0,1] ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Values in dimensionsmask should be 0 (inactive) or 1 (active)."
                    f" At position {str(idm)} given {str(dm)}."
                )
        self.dimensionsmask = dimensionsmask

        # randomgenerator
        if ( not isinstance(randomgenerator, int ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for randomgenerator parameter. It should be int, but"
                f" {str(type(randomgenerator))} was given."
            )
        if ( randomgenerator not in [0,1,2] ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid value for randomgenerator. Allowed values are 0 (determined by compiler)"
                f" 1 (based on random_number) or 2 (Ziggurat method). "
                f"{str(randomgenerator)} was given."
            )
        self.randomgenerator = randomgenerator


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

        with open(self.fn_path, 'w') as f:

            # Write time step selection method
            if self.timestep is not None:
                if self.timestep == 'adv': 
                    f.write(f"ADV\n")  
                    f.write(f"{self.courant:.10f}\n")  
                if self.timestep == 'disp': 
                    f.write(f"DISP\n")
                    f.write(f"{self.ctdisp:.10f}\n")  
                if self.timestep == 'min':
                    f.write(f"MIN_ADV_DISP\n")  
                    f.write(f"{self.courant:.10f}\n")  
                    f.write(f"{self.ctdisp:.10f}\n")  
                if self.timestep == 'fixed':
                    f.write(f"FIXED\n")  
                    f.write(f"{self.deltat:.10f}\n")  

            # Write advection model
            f.write(f"{self.advection.upper():20s}\n")  

            # dimensionsmask
            for idb, b in enumerate(self.dimensionsmask):
                if idb == len(self.dimensionsmask)-1: 
                    f.write(f"{b:10d}\n")
                else:
                    f.write(f"{b:10d}\t")

            # randomgenerator
            if ( self.randomgenerator != 0 ): 
                f.write(f"{self.randomgenerator}\n")

        # Done
        return
