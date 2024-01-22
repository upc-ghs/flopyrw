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


class dispersionModel(Enum): 
    '''
    Enumerate input formats for the package
    '''
    iso          = 0
    isotropic    = 0
    axi          = 1
    axisymmetric = 1


class ModpathRWDsp( Package ):
    """
    MODPATH-RW Dispersion Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class: ModpathRW) to which
        this package will be added.
    modelkind : str
        The dispersion model. Determines interpretation of dispersivities. 
        Allowed values are: 
          * iso/isotropic: transverse dispersion is isotropic, interpret parameters
            alphal, alphat. 
          * axi/axisymmetric: the dispersion model from Lichtner et al. 2002 for 
            vertical symmetry axis. Interpret four dispersivities: alphalh, alphalv,
            alphath and alphatv. Index h means perpendicular to the axis of symmetry 
            and index v means parallel to the axis of symmetry.
    inputformat : str
        The format for writing dispersion model paramters. Allowed values are: 
          * uni/uniform: Parameters are spatially uniform
          * dist/distributed/u3d: Parameters are spatially distributed
        Notes:
         * The inputformat distinction determines the reader to be used in MODPATH-RW.
           In the case of spatially uniform parameters, arrays are allocated with only 
           one value instead of the default (one value per cell) while using u3d readers.
         * If None is given, it will infer an inputformat from the shape of dispersion parameters.
    alphal : float, np.ndarray
        Longitudinal dispersivity. Should satisfy alphal >=0.
        Interpreted when modelkind = iso/isotropic
    alphat : float, np.ndarray
        Transverse dispersivity. Should satisfy alphat >=0.
        Interpreted when modelkind = iso/isotropic
    alphalh : float, np.ndarray
        Longitudinal dispersivity for horizontal flow. Should satisfy alphalh >=0.
        Interpreted when modelkind = axi/axisymmetric
    alphalv : float, np.ndarray
        Longitudinal dispersivity for vertical flow. Should satisfy alphalv >=0.
        Interpreted when modelkind = axi/axisymmetric
    alphath : float, np.ndarray
        Horizontal transverse dispersivity. Should satisfy alphath >=0.
        Interpreted when modelkind = axi/axisymmetric
    alphatv : float, np.ndarray
        Vertical transverse dispersivity. Should satisfy alphatv >=0.
        Interpreted when modelkind = axi/axisymmetric
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
        modelkind   = 'iso'   , # iso/isotropic, axi/axisymmetric
        inputformat = None    , # uni/uniform or dist/distributed/u3d
        alphal      = 0.1     , # xx
        alphat      = 0.01    , # yy, zz isotropic transverse dispersivity
        dmeff       = 0.0     , # corrected by tortuosity
        alphalh     = 0.1     , # longitudinal for horizontal flow 
        alphalv     = 0.1     , # longitudinal for vertical flow
        alphath     = 0.01    , # yy, zz isotropic transverse dispersivity
        alphatv     = 0.01    , # yy, zz isotropic transverse dispersivity
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
                model.multipackage[ftype]['unitnumber'],
                allowDuplicates=True,
            )
        else:
            # Pass the parent
            self._parent = model.multipackage[ftype]['instances'][0]._parent


        # modelkind
        if ( not isinstance( modelkind, str ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for modelkind. It should be a str , but"
                f" {str(type(modelkind))} was given."
            )
        if modelkind is not None: 
            try:
                self.modelkindid = dispersionModel[modelkind.lower()].value
                self.modelkind   = modelkind.lower()
            except:
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid dispersion modelkind {str(modelkind)}."
                    f" The allowed values are iso/isotropic or axi/axisymmetric."
                )
        else:
            self.modelkind   = 'iso'
            self.modelkindid = 0


        # basic preliminary health checks for dispersion params params
        if self.modelkindid == 0:
            if ( not isinstance(alphal,(int,float,np.ndarray,list)) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for alphal. It should be a real value int/float/list/np.ndarray."
                    f" {str(type(alphal))} was given."
                )
            if ( not isinstance(alphat,(int,float,np.ndarray,list)) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for alphat. It should be a real value int/float/list/np.ndarray."
                    f" {str(type(alphat))} was given."
                )
            if ( not isinstance(dmeff,(int,float,np.ndarray,list)) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for dmeff. It should be a real value int/float/list/np.ndarray."
                    f" {str(type(dmeff))} was given."
                )
        elif self.modelkindid == 1:
            if ( not isinstance(alphalh,(int,float,np.ndarray,list)) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for alphalh. It should be a real value int/float/list/np.ndarray."
                    f" {str(type(alphalh))} was given."
                )
            if ( not isinstance(alphalv,(int,float,np.ndarray,list)) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for alphalv. It should be a real value int/float/list/np.ndarray."
                    f" {str(type(alphalv))} was given."
                )
            if ( not isinstance(alphath,(int,float,np.ndarray,list)) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for alphath. It should be a real value int/float/list/np.ndarray."
                    f" {str(type(alphath))} was given."
                )
            if ( not isinstance(alphatv,(int,float,np.ndarray,list)) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for alphatv. It should be a real value int/float/list/np.ndarray."
                    f" {str(type(alphatv))} was given."
                )
            if ( not isinstance(dmeff,(int,float,np.ndarray,list)) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for dmeff. It should be a real value int/float/list/np.ndarray."
                    f" {str(type(dmeff))} was given."
                )


        # If no specific input format was given, 
        # infer from the type of parameters
        if ( inputformat is None ):
            # linear isotropic
            if self.modelkindid == 0:
                # uniform params
                if ( 
                    ( isinstance(alphal,(int,float)) ) and 
                    ( isinstance(alphat,(int,float)) ) and 
                    ( isinstance(dmeff ,(int,float)) )  
                ):
                    inputformat = 'uni'
                else:
                # non uniform params
                    if ( isinstance( alphal, (int,float) ) ):
                        alphal = [alphal]
                    if ( isinstance( alphat, (int,float) ) ):
                        alphat = [alphat]
                    if ( isinstance( dmeff, (int,float) ) ):
                        dmeff  = [dmeff]
                    if ( isinstance(alphal,list) ): 
                        alphal = np.array(alphal)
                    if ( isinstance(alphat,list) ): 
                        alphat = np.array(alphat)
                    if ( isinstance(dmeff,list) ): 
                        dmeff = np.array(dmeff)
                    if ( alphal.dtype not in [np.int32, np.int64, np.float32, np.float64] ):
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for alphal. It should be a real value, but"
                            f" {str(alphal.dtype)} was given."
                        )
                    if ( alphat.dtype not in [np.int32, np.int64, np.float32, np.float64] ):
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for alphat. It should be a real value, but"
                            f" {str(alphat.dtype)} was given."
                        )
                    if ( dmeff.dtype not in [np.int32, np.int64, np.float32, np.float64] ):
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for dmeff. It should be a real value, but"
                            f" {str(dmeff.dtype)} was given."
                        )
                    inputformat = 'dist'
            # linear axisymmetric
            elif self.modelkindid == 1:
                # uniform params
                if ( 
                    ( isinstance(alphalh,(int,float)) ) and 
                    ( isinstance(alphalv,(int,float)) ) and 
                    ( isinstance(alphath,(int,float)) ) and 
                    ( isinstance(alphatv,(int,float)) ) and 
                    ( isinstance(dmeff ,(int,float)) )  
                ):
                    inputformat = 'uni'
                else:
                # non uniform params
                    if ( isinstance( alphalh, (int,float) ) ):
                        alphalh = [alphalh]
                    if ( isinstance( alphalv, (int,float) ) ):
                        alphalv = [alphalv]
                    if ( isinstance( alphath, (int,float) ) ):
                        alphath = [alphath]
                    if ( isinstance( alphatv, (int,float) ) ):
                        alphatv = [alphatv]
                    if ( isinstance( dmeff, (int,float) ) ):
                        dmeff  = [dmeff]
                    if ( isinstance(alphalh,list) ): 
                        alphalh = np.array(alphalh)
                    if ( isinstance(alphalv,list) ): 
                        alphalv = np.array(alphalv)
                    if ( isinstance(alphath,list) ): 
                        alphath = np.array(alphath)
                    if ( isinstance(alphatv,list) ): 
                        alphatv = np.array(alphatv)
                    if ( isinstance(dmeff,list) ): 
                        dmeff = np.array(dmeff)
                    if ( alphalh.dtype not in [np.int32, np.int64, np.float32, np.float64] ):
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for alphalh. It should be a real value, but"
                            f" {str(alphalh.dtype)} was given."
                        )
                    if ( alphalv.dtype not in [np.int32, np.int64, np.float32, np.float64] ):
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for alphalv. It should be a real value, but"
                            f" {str(alphalv.dtype)} was given."
                        )
                    if ( alphath.dtype not in [np.int32, np.int64, np.float32, np.float64] ):
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for alphath. It should be a real value, but"
                            f" {str(alphath.dtype)} was given."
                        )
                    if ( alphatv.dtype not in [np.int32, np.int64, np.float32, np.float64] ):
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for alphatv. It should be a real value, but"
                            f" {str(alphatv.dtype)} was given."
                        )
                    if ( dmeff.dtype not in [np.int32, np.int64, np.float32, np.float64] ):
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for dmeff. It should be a real value, but"
                            f" {str(dmeff.dtype)} was given."
                        )
                    inputformat = 'dist'

        # inputformat
        if ( not isinstance( inputformat, str ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for inputformat option {str(type(inputformat))}."
                f" It was expecting a string." 
            )
        try:
            self.inputformat = inputFormat[inputformat.lower()].value
        except:
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid input format option {str(inputformat)}."
                f" The allowed values are uni/uniform or dis/distributed/u3d." 
            )

        # Assign params
        if self.inputformat == 0:
            # Spatially constant, uniform 
            if self.modelkindid == 0:
                # Validate types
                # alphal
                if ( not isinstance(alphal,(int,float)) ):
                    if ( isinstance(alphal,(list,np.ndarray)) and (len(alphal)==1)):
                        alphal = alphal[0]
                    else:
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for alphal. Uniform format was requested, but type"
                            f" {str(type(alphal))} was given."
                        )
                # alphat
                if ( not isinstance(alphat,(int,float)) ):
                    if ( isinstance(alphat,(list,np.ndarray)) and (len(alphat)==1)):
                        alphat = alphat[0]
                    else:
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for alphat. Uniform format was requested, but type"
                            f" {str(type(alphat))} was given."
                        )
                # dmeff
                if ( not isinstance(dmeff,(int,float)) ):
                    if ( isinstance(dmeff,(list,np.ndarray)) and (len(dmeff)==1)):
                        dmeff = dmeff[0]
                    else:
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for dmeff. Uniform format was requested, but type"
                            f" {str(type(dmeff))} was given."
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
                        f"{self.__class__.__name__}:"
                        f" Inconsistent value of alphat {str(alphat)}."
                        f" It should be greater or equal than zero."
                    )
                if ( dmeff < 0 ):
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Inconsistent value of dmeff {str(dmeff)}."
                        f" It should be greater or equal than zero."
                    )
                self.alphal = alphal 
                self.alphat = alphat
                self.dmeff  = dmeff

            elif self.modelkindid == 1:
                # linear axisymmetric
                # alphalh
                if ( not isinstance(alphalh,(int,float)) ):
                    if ( isinstance(alphalh,(list,np.ndarray)) and (len(alphalh)==1)):
                        alphalh = alphalh[0]
                    else:
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for alphalh. Uniform format was requested, but type"
                            f" {str(type(alphalh))} was given."
                        )
                # alphalv
                if ( not isinstance(alphalv,(int,float)) ):
                    if ( isinstance(alphalv,(list,np.ndarray)) and (len(alphalv)==1)):
                        alphalv = alphalv[0]
                    else:
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for alphalv. Uniform format was requested, but type"
                            f" {str(type(alphalv))} was given."
                        )
                # alphath
                if ( not isinstance(alphath,(int,float)) ):
                    if ( isinstance(alphath,(list,np.ndarray)) and (len(alphath)==1)):
                        alphath = alphath[0]
                    else:
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for alphath. Uniform format was requested, but type"
                            f" {str(type(alphath))} was given."
                        )
                # alphatv
                if ( not isinstance(alphatv,(int,float)) ):
                    if ( isinstance(alphatv,(list,np.ndarray)) and (len(alphatv)==1)):
                        alphatv = alphatv[0]
                    else:
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for alphatv. Uniform format was requested, but type"
                            f" {str(type(alphatv))} was given."
                        )
                # dmeff
                if ( not isinstance(dmeff,(int,float)) ):
                    if ( isinstance(dmeff,(list,np.ndarray)) and (len(dmeff)==1)):
                        dmeff = dmeff[0]
                    else:
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid type for dmeff. Uniform format was requested, but type"
                            f" {str(type(dmeff))} was given."
                        )
                # Validate values
                if ( alphalh < 0 ):
                    raise ValueError(
                        self.__class__.__name__ + ':'
                        + ' Inconsistent value of alphalh ' +str(alphalh)
                        + '. It should be greater or equal than zero.'
                    )
                if ( alphalv < 0 ):
                    raise ValueError(
                        self.__class__.__name__ + ':'
                        + ' Inconsistent value of alphalv ' +str(alphalv)
                        + '. It should be greater or equal than zero.'
                    )
                if ( alphath < 0 ):
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Inconsistent value of alphath {str(alphath)}."
                        f" It should be greater or equal than zero."
                    )
                if ( alphatv < 0 ):
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Inconsistent value of alphatv {str(alphatv)}."
                        f" It should be greater or equal than zero."
                    )
                if ( dmeff < 0 ):
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Inconsistent value of dmeff {str(dmeff)}."
                        f" It should be greater or equal than zero."
                    )
                self.alphalh = alphalh
                self.alphalv = alphalv 
                self.alphath = alphath
                self.alphatv = alphatv
                self.dmeff   = dmeff
            else:
                raise NotImplementedError(
                        f"{self.__class__.__name__}:"
                        f" Dispersion model kind {str(self.modelkind)}"
                        f" is not implemented."
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

            if self.modelkindid == 0:
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
                        f"{self.__class__.__name__}:"
                        f" Error while initializing distributed variable alphal."
                        f" Is the input shape consistent with flow model dimensions ?"
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
                        f"{self.__class__.__name__}:"
                        f" Error while initializing distributed variable alphat."
                        f" Is the input shape consistent with flow model dimensions ?"
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
                        f"{self.__class__.__name__}:"
                        f" Error while initializing distributed variable dmeff."
                        f" Is the input shape consistent with flow model dimensions ?"
                    )

                # Validate values
                if np.any( self.alphal.array < 0 ) : 
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Inconsistent values for alphal. Negatives were found and"
                        f" it should be greater or equal than zero."
                    )
                if np.any( self.alphat.array < 0 ) : 
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Inconsistent values for alphat. Negatives were found and"
                        f" it should be greater or equal than zero."
                    )
                if np.any( self.dmeff.array < 0 ) : 
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Inconsistent values for dmeff. Negatives were found and"
                        f" it should be greater or equal than zero."
                    )
            elif self.modelkindid == 1:
                # linear axisymmetric
                # alphalh
                try:
                    self.alphalh = Util3d(
                        model,
                        shape3d,
                        np.float32,
                        alphalh,
                        name="ALPHALH",
                        locat=self._parent.multipackage[ftype]['unitnumber'], 
                    )
                except:
                    raise Exception(
                        f"{self.__class__.__name__}:"
                        f" Error while initializing distributed variable alphalh."
                        f" Is the input shape consistent with flow model dimensions ?"
                    )
                # alphalv
                try:
                    self.alphalv = Util3d(
                        model,
                        shape3d,
                        np.float32,
                        alphalv,
                        name="ALPHALV",
                        locat=self._parent.multipackage[ftype]['unitnumber'], 
                    )
                except:
                    raise Exception(
                        f"{self.__class__.__name__}:"
                        f" Error while initializing distributed variable alphalv."
                        f" Is the input shape consistent with flow model dimensions ?"
                    )
                # alphath
                try:
                    self.alphath = Util3d(
                        model,
                        shape3d,
                        np.float32,
                        alphath,
                        name="ALPHATH",
                        locat=self._parent.multipackage[ftype]['unitnumber'], 
                    )
                except:
                    raise Exception(
                        f"{self.__class__.__name__}:"
                        f" Error while initializing distributed variable alphath."
                        f" Is the input shape consistent with flow model dimensions ?"
                    )
                # alphatv
                try:
                    self.alphatv = Util3d(
                        model,
                        shape3d,
                        np.float32,
                        alphatv,
                        name="ALPHATV",
                        locat=self._parent.multipackage[ftype]['unitnumber'], 
                    )
                except:
                    raise Exception(
                        f"{self.__class__.__name__}:"
                        f" Error while initializing distributed variable alphatv."
                        f" Is the input shape consistent with flow model dimensions ?"
                    )
                # dmeff
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
                        f"{self.__class__.__name__}:"
                        f" Error while initializing distributed variable dmeff."
                        f" Is the input shape consistent with flow model dimensions ?"
                    )

                # Validate values
                if np.any( self.alphalh.array < 0 ) : 
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Inconsistent values for alphalh. Negatives were found and"
                        f" it should be greater or equal than zero."
                    )
                if np.any( self.alphalv.array < 0 ) : 
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Inconsistent values for alphalv. Negatives were found and"
                        f" it should be greater or equal than zero."
                    )
                if np.any( self.alphath.array < 0 ) : 
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Inconsistent values for alphath. Negatives were found and"
                        f" it should be greater or equal than zero."
                    )
                if np.any( self.alphatv.array < 0 ) : 
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Inconsistent values for alphatv. Negatives were found and"
                        f" it should be greater or equal than zero."
                    )
                if np.any( self.dmeff.array < 0 ) : 
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Inconsistent values for dmeff. Negatives were found and"
                        f" it should be greater or equal than zero."
                    )
            else:
                raise NotImplementedError(
                        f"{self.__class__.__name__}:"
                        f" Dispersion model kind {str(self.modelkind)}"
                        f" is not implemented."
                    )

        # String id 
        if (stringid is not None): 
            if ( not isinstance( stringid, str ) ): 
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for stringid. It should be str, but"
                    f" {str(type(stringid))} was given."
                )
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
                if ins.modelkindid == 0:
                # linear isotropic
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
                if ins.modelkindid == 1:
                # linear axisymmetric
                    if ins.inputformat == 0:
                        # Uniform
                        f.write(f"{ins.dmeff :.10f}\n")  
                        f.write(f"{ins.alphalh:.10f}\n")  
                        f.write(f"{ins.alphalv:.10f}\n")  
                        f.write(f"{ins.alphath:.10f}\n")  
                        f.write(f"{ins.alphatv:.10f}\n")  
                    elif ins.inputformat == 1:
                        # Distributed
                        f.write(ins.dmeff.get_file_entry())
                        f.write(ins.alphalh.get_file_entry())
                        f.write(ins.alphalv.get_file_entry())
                        f.write(ins.alphath.get_file_entry())
                        f.write(ins.alphatv.get_file_entry())
       
        # Done 
        return
