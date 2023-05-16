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
            'adv'   : local time step computed with advection criteria
            'disp'  : local time step computed with dispersion criteria
            'min'   : local time step obtained as the minimum between 'adv' and 'disp'
            'fixed' : the same for all particles, given by the user in the 'deltat' keyword 
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
    extension      : str, optional
        File extension (default is 'rwopts').
    """

    def __init__(
        self,
        model,
        timestep         = 'min',
        courant          = 0.1,
        ctdisp           = 0.1,
        deltat           = 1.0,
        advection        = 'eulerian',
        dimensionsmask   = [1,1,1],
        extension        = 'rwopts',
    ):

        unitnumber = model.next_unit()

        super().__init__(model, extension, "RWOPTS", unitnumber)

        # timestep
        if ( timestep not in ['adv', 'disp', 'min', 'fixed'] ):
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Time step selection ' +str(timestep)+ ' is not valid.'+
                '. Allowed values are: adv, disp, min, fixed.'
            )
        self.timestep = timestep

        # timestep selection params
        if (courant>0): 
            self.courant  = courant
        else:
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Courant number for advection criteria should be greater than zero.'
            )
        if (ctdisp>0): 
            self.ctdisp  = ctdisp
        else:
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' CT constant for dispersion criteria should be greater than zero.'
            )
        if (deltat>0): 
            self.deltat  = deltat
        else:
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Time step for particles displacement should be greater than zero.'
            )

        # advection 
        if ( advection not in ['eulerian','exponential'] ):
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Advection model ' +str(timestep)+ ' is not valid.'+
                '. Allowed values are: eulerian, exponential.'
            )
        self.advection  = advection

        # dimensionsmask
        if not isinstance( dimensionsmask, (list, np.array) ):
            raise TypeError(
                self.__class__.__name__ + ':' + 
                ' dimensionsmask should be of type list or np.array.' +
                '. Given ' + str(type(dimensionsmask))
            )
        if isinstance( dimensionsmask, list ):
            dimensionsmask = np.array(dimensionsmask).astype(np.int32)
        if len(dimensionsmask) != 3:
            if len(dimensionsmask) != 3:
                raise ValueError(
                    self.__class__.__name__ + ':' + 
                    ' dimensionsmask should be of length 3.' +
                    '. Given ' + str(len(dimensionsmask))
                )
        for nd in range(3):
            if (dimensionsmask[nd] not in [0,1] ):
                raise ValueError(
                    self.__class__.__name__ + ':' + 
                    ' values in dimensionsmask should be 0 (inactive) or 1 (active).' +
                    ' At position ', str(nd), ' given ', str(dimensionsmask[nd])
                )
        self.dimensionsmask = dimensionsmask


        self.parent.add_package(self)


        # Done !
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
        f = open(self.fn_path, "w")

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

        # And close
        f.close()

        return
        
