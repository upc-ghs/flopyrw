'''
Configuration of MODPATH-RW dispersion
'''

# python
import numpy as np
from enum import Enum

# flopy
from flopy.pakbase import Package
from flopy.utils import Util3d

# local 
from .utils import multipackage

class inputFormat(Enum): 
    '''
    Enumerate input formats for the package
    '''
    uni         = 0
    uniform     = 0
    dist        = 1
    distributed = 1
    u3d         = 1


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
    inputformat : str
        The format for writing dispersion model paramters. Allowed values are: 
            uni/uniform: Parameters are spatially uniform
            dist/distributed/u3d: Parameters are spatially distributed
        Note: this distinction determines the reader to be used in MODPATH-RW. In the case of 
        spatially uniform parameters, arrays are allocated with only one value instead of the 
        default (one value per cell) while using u3d readers.
    alphal : float, np.ndarray
        Longitudinal dispersivity. Should satisfy alphal >=0.
    alphat : float, np.ndarray
        Transverse dispersivity. Should satisfy alphat >=0.
    dmeff  : float, np.array
        Effective molecular diffusion, corrected by tortuosity. Should satisfy dmeff >=0.
    stringid : str, optional
        An id for the dispersion specification. By default concatenates
        package ftype and integer id. 
    extension : str, optional
        File extension (default is 'dsp').
    """


    @staticmethod
    def _ftype():
        return 'DSP'


    @multipackage
    def __init__(
        self,
        model,
        modelkind   = 'linear', # linear
        inputformat = None    , # uni/uniform or dist/distributed/u3d
        alphal      = 0.1     , # xx
        alphat      = 0.01    , # yy, zz isotropic transverse dispersivity
        dmeff       = 0.0     , # corrected by tortuosity
        stringid    = None    ,
        extension   = 'dsp'   ,
    ):

        # Call parent constructor
        ftype = self._ftype()
        if ( model.multipackage[ftype]['count'] == 0 ): 
            super().__init__(
                model,
                extension,
                ftype,
                model.multipackage[ftype]['unitnumber']
            )
        else:
            # Pass the parent
            self._parent = model.multipackage[ftype]['instances'][0]._parent


        # modelkind
        if modelkind is not None: 
            if modelkind not in ('linear'):
                raise ValueError(
                    self.__class__.__name__ + ':'
                    + ' Invalid kind option ' + str(modelkind)
                    + '. The only allowed value is linear. '
                )
            if modelkind == 'linear':
                self.modelkind   = modelkind
                self.modelkindid = 1
        else:
            self.modelkind   = 'linear'
            self.modelkindid = 1


        # basic preliminary health checks for dispersion params params
        if self.modelkindid == 1:
            if ( not isinstance(alphal,(int,float,np.ndarray,list)) ):
                raise TypeError(
                    self.__class__.__name__ + ':'
                    + ' Invalid type for alphal. It should be a real value int/float/list/np.ndarray. '
                    + str(type(alphal)) + ' was given. '
                )
            if ( not isinstance(alphat,(int,float,np.ndarray,list)) ):
                raise TypeError(
                    self.__class__.__name__ + ':'
                    + ' Invalid type for alphat. It should be a real value int/float/list/np.ndarray. '
                    + str(type(alphat)) + ' was given. '
                )
            if ( not isinstance(alphat,(int,float,np.ndarray,list)) ):
                raise TypeError(
                    self.__class__.__name__ + ':'
                    + ' Invalid type for dmeff. It should be a real value int/float/list/np.ndarray. '
                    + str(type(dmeff)) + ' was given. '
                )


        # If no specific input format was given, 
        # infer from the type of parameters
        if ( inputformat is None ):
            # uniform params
            if ( 
                ( isinstance(alphal,(int,float)) ) and 
                ( isinstance(alphat,(int,float)) ) and 
                ( isinstance(dmeff ,(int,float)) )  
            ):
                inputformat = 'uni'
            else:
            # non uniform params
                if ( isinstance(alphal,list) ): 
                    alphal = np.array(alphal)
                if ( isinstance(alphat,list) ): 
                    alphat = np.array(alphat)
                if ( isinstance(dmeff,list) ): 
                    dmeff = np.array(dmeff)
                if ( alphal.dtype not in [np.int32, np.int64, np.float32, np.float64] ):
                    raise TypeError(
                        self.__class__.__name__ + ':' 
                        + ' Invalid type for alphal. It should be a real value, but '
                        + str(alphal.dtype) + ' was given. '
                    )
                if ( alphat.dtype not in [np.int32, np.int64, np.float32, np.float64] ):
                    raise TypeError(
                        self.__class__.__name__ + ':'
                        + ' Invalid type for alphat. It should be a real value, but '
                        + str(alphat.dtype) + ' was given. '
                    )
                if ( dmeff.dtype not in [np.int32, np.int64, np.float32, np.float64] ):
                    raise TypeError(
                        self.__class__.__name__ + ':'
                        + ' Invalid type for dmeff. It should be a real value, but '
                        + str(dmeff.dtype) + ' was given. '
                    )
                inputformat = 'dist'

        # inputformat
        try:
            self.inputformat = inputFormat[inputformat.lower()].value
        except:
            raise ValueError(
                self.__class__.__name__ + ':'
                + ' Invalid input format option ' +str(inputformat)
                + '. The allowed values are uni/uniform or dis/distributed/u3d.' 
            )

        # Assign params
        if self.inputformat == 0:
            # Spatially constant, uniform 
            if self.modelkindid == 1:
                # Validate types
                if ( not isinstance(alphal,(int,float)) ):
                    if ( isinstance(alphal,(list,np.ndarray)) and (len(alphal)==1)):
                        alphal = alphal[0]
                    else:
                        raise TypeError(
                            self.__class__.__name__ + ':'
                            + ' Invalid type for alphal. Uniform format was requested but type '
                            + str(type(alphal)) + ' was given. '
                        )
                if ( not isinstance(alphat,(int,float)) ):
                    if ( isinstance(alphat,(list,np.ndarray)) and (len(alphat)==1)):
                        alphat = alphat[0]
                    else:
                        raise TypeError(
                            self.__class__.__name__ + ':'
                            + ' Invalid type for alphat. Uniform format was requested but type '
                            + str(type(alphat)) + ' was given. '
                        )
                if ( not isinstance(dmeff,(int,float)) ):
                    if ( isinstance(dmeff,(list,np.ndarray)) and (len(dmeff)==1)):
                        dmeff = dmeff[0]
                    else:
                        raise TypeError(
                            self.__class__.__name__ + ':'
                            + ' Invalid type for dmeff. Uniform format was requested but type '
                            + str(type(dmeff)) + ' was given. '
                        )
                # Validate values
                if ( alphal < 0 ):
                    raise ValueError(
                        self.__class__.__name__ + ':'
                        + ' Inconsistent value of alphal ' +str(alphal)
                        + '. It should be greater or equal than zero.'
                    )
                if ( alphat < 0 ):
                    raise ValueError(
                        self.__class__.__name__ + ':'
                        + ' Inconsistent value of alphat ' +str(alphat)
                        + '. It should be greater or equal than zero.'
                    )
                if ( dmeff < 0 ):
                    raise ValueError(
                        self.__class__.__name__ + ':'
                        + ' Inconsistent value of dmeff ' +str(dmeff)
                        + '. It should be greater or equal than zero.'
                    )
                self.alphal = alphal 
                self.alphat = alphat
                self.dmeff  = dmeff
            else:
                raise NotImplementedError(
                        self.__class__.__name__ + ':'
                        + ' Dispersion model kind ' +str(self.modelkind)
                        + ' is not implemented.'
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
                try:
                    self.alphal = Util3d(
                        model,
                        shape3d,
                        np.float32,
                        alphal,
                        name="ALPHAL",
                        locat=self._parent.multipackage[ftype]['unitnumber'], 
                    )
                except:
                    raise Exception(
                        self.__class__.__name__ + ':'
                        + ' Error while initializing distributed variable alphal. '
                        + 'Is the input shape consistent with flow model dimensions ? '
                    )
                try:
                    self.alphat = Util3d(
                        model,
                        shape3d,
                        np.float32,
                        alphat,
                        name="ALPHAT",
                        locat=self._parent.multipackage[ftype]['unitnumber'], 
                    )
                except:
                    raise Exception(
                        self.__class__.__name__ + ':'
                        + ' Error while initializing distributed variable alphat. '
                        + 'Is the input shape consistent with flow model dimensions ? '
                    )
                try:
                    self.dmeff = Util3d(
                        model,
                        shape3d,
                        np.float32,
                        dmeff,
                        name="DMEFF",
                        locat=self._parent.multipackage[ftype]['unitnumber'], 
                    )
                except:
                    raise Exception(
                        self.__class__.__name__ + ':'
                        + ' Error while initializing distributed variable dmeff. '
                        + 'Is the input shape consistent with flow model dimensions ? '
                    )

                # Validate values
                if np.any( self.alphal.array < 0 ) : 
                    raise ValueError(
                        self.__class__.__name__ + ':'
                        + ' Inconsistent values for alphal. Negatives were found and '
                        + 'it should be greater or equal than zero.'
                    )
                if np.any( self.alphat.array < 0 ) : 
                    raise ValueError(
                        self.__class__.__name__ + ':'
                        + ' Inconsistent values for alphat. Negatives were found and '
                        + 'it should be greater or equal than zero.'
                    )
                if np.any( self.dmeff.array < 0 ) : 
                    raise ValueError(
                        self.__class__.__name__ + ':'
                        + ' Inconsistent values for dmeff. Negatives were found and '
                        + 'it should be greater or equal than zero.'
                    )
            else:
                raise NotImplementedError(
                        self.__class__.__name__ + ':'
                        + ' Dispersion model kind ' +str(self.modelkind)
                        + ' is not implemented.'
                    )

        # String id 
        if (stringid is not None): 
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
        return
