'''
Configuration of MODPATH-RW dispersion
'''

# python
import numpy as np

# flopy
from flopy.pakbase import Package
from flopy.utils import Util3d

# local 
from .utils import count_instances # Increment COUNTER


class ModpathRWDsp( Package ):
    """
    MODPATH-RW Dispersion Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class: ModpathRW) to which
        this package will be added.
    modelkind : str
        The dispersion model. The only currently implemented is with linear dispersivities 
        and isotropic transverse dispersivity. Allowed values 'linear'.
    inputformat : int
        The format for writing dispersion model paramters. Allowed values are: 
            0: Parameters are spatially uniform
            1: Parameters are spatially distributed
        Note: this distinction determines the reader to be used in MODPATH-RW. In the case of 
        spatially uniform parameters, arrays are allocated with only one value instead of the 
        default (one value per cell) while using u3d readers.
    alphal : float
        Longitudinal dispersivity. Should satisfy alphal >=0.
    alphat : float
        Transverse dispersivity. Should satisfy alphat >=0.
    dmeff  : float
        Effective molecular diffusion, corrected by tortuosity. Should satisfy dmeff >=0.
    id : int 
        A positive integer identifier. If None given is automatically assigned with 
        the instance counter. 
    stringid : str, optional
        An id for the dispersion specification. If not given is automatically filled with 
        DSP+str(id), where id is an internal use counter. 
    extension : str, optional
        File extension (default is 'dsp').
    """

    # Class properties
    COUNTER    = 0
    UNITNUMBER = 0
    INSTANCES  = []

    @count_instances
    def __init__(
        self,
        model,
        modelkind         = 'linear', # linear
        inputformat       = 0       , # 0: uniform, 1: u3d
        alphal            = 0.1     , # xx
        alphat            = 0.01    , # yy, zz isotropic transverse dispersivity
        dmeff             = 0.0     , # corrected by tortuosity
        alphath           = 0.01    , # ( not implemented ) 
        alphatv           = 0.01    , # ( not implemented )
        dmaqueous         = 0.0     , # ( not implemented )  
        betal             = 1       , # ( not implemented ) 
        betath            = 0.5     , # ( not implemented )  
        betatv            = 0.5     , # ( not implemented )  
        delta             = 5       , # ( not implemented )  
        dgrain            = 1       , # ( not implemented ) 
        id                = None    , 
        stringid          = None    ,
        extension         = 'dsp'   ,
    ):

        # Define UNITNUMBER if the first instance is created
        if self.__class__.COUNTER == 1:
            unitnumber = model.next_unit()
            super().__init__(model, extension, "DSP", unitnumber)
            self.__class__.UNITNUMBER = self.unit_number[0]
        else:
            # Needed for Util3d
            self._parent = self.INSTANCES[0]._parent

        # Select modelkind
        if modelkind is not None: 
            if modelkind not in ('linear'):
                raise ValueError(
                    self.__class__.__name__ + ':' + 
                    ' Invalid kind option ' +str(modelkind)+
                    '. The only allowed value is linear. '
                )
            if modelkind == 'linear':
                self.modelkind   = modelkind
                self.modelkindid = 1
            elif modelkind == 'nonlinear':
                self.modelkind   = modelkind
                self.modelkindid = 2
        else:
            self.modelkind   = 'linear'
            self.modelkindid = 1

        if (inputformat not in [0,1]):
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid input format option ' +str(inputformat)+
                '. The allowed values are 0 (uniform) or 1 (distributed, u3d).'
            )
        self.inputformat = inputformat 

        # Assign params
        if self.inputformat == 0:
            # Needs some handling for the case of array input,
            # extract the first value or something.

            # Spatially constant, uniform 
            if self.modelkindid == 1:
                if ( ( alphal is None ) or (alphal < 0) ):
                    raise ValueError(
                        self.__class__.__name__ + ':' + 
                        ' Inconsitent value of alphal ' +str(alphal)+
                        '. It should be greater or equal than zero.'
                    )
                if ( ( alphat is None ) or (alphat < 0) ):
                    raise ValueError(
                        self.__class__.__name__ + ':' + 
                        ' Inconsitent value of alphat ' +str(alphat)+
                        '. It should be greater or equal than zero.'
                    )
                if ( ( dmeff is None ) or (dmeff < 0) ):
                    raise ValueError(
                        self.__class__.__name__ + ':' + 
                        ' Inconsitent value of dmeff ' +str(dmeff)+
                        '. It should be greater or equal than zero.'
                    )
                self.alphal = alphal 
                self.alphat = alphat
                self.dmeff  = dmeff
            else:
                raise NotImplementedError(
                        self.__class__.__name__ + ':' + 
                        ' Dispersion model kind ' +str(self.modelkind)+
                        '. Is not implemented.'
                    )
        elif self.inputformat == 1:
            # Distributed dispersion parameters 
            shape = model.shape
            if len(shape) == 3:
                shape3d = shape
            elif len(shape) == 2:
                shape3d = (shape[0], 1, shape[1])
            else:
                shape3d = (1, 1, shape[0])
            self.model = model
            self.shape3d = shape3d

            if self.modelkindid == 1:
                # linear isotropic
                self.alphal = Util3d(
                    model,
                    shape3d,
                    np.float32,
                    alphal,
                    name="ALPHAL",
                    locat=self.__class__.UNITNUMBER,
                )
                self.alphat = Util3d(
                    model,
                    shape3d,
                    np.float32,
                    alphath,
                    name="ALPHAT",
                    locat=self.__class__.UNITNUMBER,
                )
                self.dmeff = Util3d(
                    model,
                    shape3d,
                    np.float32,
                    dmeff,
                    name="DMEFF",
                    locat=self.__class__.UNITNUMBER,
                )
            else:
                raise NotImplementedError(
                        self.__class__.__name__ + ':' + 
                        ' Dispersion model kind ' +str(self.modelkind)+
                        '. Is not implemented.'
                    )

        # Define id
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER

        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = 'DSP'+str(self.__class__.COUNTER)

        # Add package and save instance        
        if self.__class__.COUNTER == 1:
            self.parent.add_package(self)
        self.__class__.INSTANCES.append( self )

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
            f.write(f"{ins.stringid}\n")

            # Write modelkindid and input format
            f.write(f"{ins.modelkindid}    {ins.inputformat}\n")
           
            # Write dispersion model parameters
            if ins.modelkindid == 1:
                if ins.inputformat == 0:
                    # Uniform
                    f.write(f"{ins.dmeff :.10f}\n")  
                    f.write(f"{ins.alphal:.10f}\n")  
                    f.write(f"{ins.alphat:.10f}\n")  
                elif ins.inputformat == 1:
                    # Distributed
                    f.write(ins.dmeff.get_file_entry())
                    f.write(ins.alphal.get_file_entry())
                    f.write(ins.alphat.get_file_entry())
       
        # Done
        f.close()


        return
