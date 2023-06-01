'''
Configuration of MODPATH-RW initial conditions
'''

# python
import numpy as np

# flopy
from flopy.pakbase import Package
from flopy.utils import Util3d

# local 
from .utils import multipackage


class ModpathRWIc( Package ):
    """
    MODPATH-RW Initial Condition Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class: ModpathRW) to which
        this package will be added.
    kind  : int 
        The kind of initial condition specification. All formats load concentration 
        with u3d. Allowed values are: 
          * 0: concentration is transformed into a distribution of particles, 
               by considering the cell volume, porosity and retardation factor.
               Particles mass is recomputed in order to impose consistency between 
               the total mass contained in the distribution and the estimated with 
               particles.
          * 1: similar to the previous alternative, particlesmass is employed to estimate 
               a number of particles per cell. The final mass of particles is obtained 
               on a per cell basis, by imposing consistency of total mass locally.
          * 2: all cells contain the same number of particles. Particles mass is calculated
               from the total solute mass at each cell.  
    particlesmass  : float
        The scale of particles mass used for transforming the concentration data into 
        a distribution of particles. 
    concentration : list, ndarray, int or float
        Interpreted as a resident initial concentration. It will be written with Util3d, 
        meaning that if the input is multidimensional, it should be consistent with 
        flow model dimensions. 
    speciesid : int 
        The species id to which this initial concentration is related. Will only 
        be written if particlesmassoption == 2. It should be given as the zero-based 
        index of the specie (the python-index of the specie in the list of species), 
        but is written to the package file with fortran format (+1, see the write_file method).
    stringid : str
        A string identifier. If None is given is automatically filled by 
        concatenating the package ftype and an integer id. 
    extension : str
        The file extension for the package. Defaults to ic. 
    """


    @staticmethod
    def _ftype():
        return 'IC'


    @multipackage
    def __init__(
        self,
        model,
        kind          = 0   ,
        particlesmass = 1.0 ,
        concentration = 0.0 ,
        speciesid     = 0   , 
        stringid      = None,
        extension     = 'ic',
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

        # Determines format for writing the initial condition
        if (kind not in [0,1]): 
        #if (kind not in [0]): 
            raise ValueError(
                self.__class__.__name__ + ':'
                + ' Invalid kind parameter ' + str(kind)
                + '. Allowed values are 0 (read concentration with u3d).'
            )
        self.kind = kind

        # The particles mass
        if ( not isinstance( particlesmass, (int,float) ) ):
            raise TypeError(
                self.__class__.__name__ + ':' 
                + ' Invalid type for particlesmass. It should be int/float, '
                + ' but ' + str(type(particlesmass)) + ' was given.' 
            )
        if (particlesmass <= 0):
            raise ValueError(
                self.__class__.__name__ + ':' 
                + ' Invalid particle mass for initial condition. '
                + 'It should be a positive value.' 
            )
        self.particlesmass = particlesmass

        # A species id, will only be written 
        # if particlesmassoption == 2
        if ( speciesid is not None ):
            if not isinstance( speciesid, int ): 
                raise TypeError(
                    self.__class__.__name__ + ':' 
                    + ' Invalid type for speciesid. It should be integer, '
                    + ' but ' + str(type(speciesid)) + ' was given.' 
                )
            if (speciesid < 0):
                raise ValueError(
                    self.__class__.__name__ + ':' 
                    + ' Invalid value for speciesid. It should be ' 
                    + ' a valid zero-based index of species.' 
                )
        self.speciesid = speciesid

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

        # and the concentration: needs id
        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])
        self.model = model
        self.shape3d = shape3d

        if (
            ( concentration is not None ) and 
            ( self.kind in [0,1] )
            #( self.kind == 0 )
        ):
            # TypeError if any other unsupported type
            if ( not isinstance( concentration, (list,float,int,np.ndarray) ) ):
                raise TypeError(
                    self.__class__.__name__ + ':'
                    + ' Invalid type ' + str(type(concentration))
                    + ' for concentration. It should be list, np.ndarray or a single numeric value.'
                )
            # If list, to np.array
            if ( isinstance( concentration, list ) ): 
                concentration = np.array( concentration )
            # Verify a valid dtype
            if ( isinstance( concentration, np.ndarray ) ): 
                if ( 
                    (concentration.dtype not in [np.int32, np.int64, np.float32, np.float64])
                ):
                    raise TypeError(
                        self.__class__.__name__ + ':' 
                        + ' Invalid dtype ' + str(concentration.dtype) 
                        + ' for concentration. It should be int or float. ' 
                    )
            try:
                self.concentration = Util3d(
                    model,
                    shape3d,
                    np.float32,
                    concentration,
                    name=self.stringid,
                    locat=self._parent.multipackage[ftype]['unitnumber'], 
                )
            except:
                raise Exception(
                    self.__class__.__name__ + ':'
                    + ' Error while initializing distributed variable concentration. '
                    + 'Is the input shape consistent with flow model dimensions ? '
                )

        elif ( (self.kind == 0) and (concentration is None) ):
            raise ValueError(
                self.__class__.__name__ + ':'
                + ' Invalid initial concentration. None was given.'
            )
    

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

                # Write string id 
                f.write(f"{ins.stringid}\n")

                # Kind/format of initial condition 
                f.write(f"{ins.kind}\n")

                # 0: resident concentration array  
                if ( ins.kind in [0,1] ):
                #if ins.kind == 0:

                    # Give particles mass 
                    f.write(f"{ins.particlesmass:.10f}\n")

                    # If solute ids are being specified, write the id
                    if( 
                        ( self._parent.particlesmassoption == 2 ) 
                    ):
                        if ( 
                            ( ins.speciesid is None ) or
                            ( ins.speciesid < 0 ) 
                        ):
                            raise ValueError(
                                self.__class__.__name__ + ':'
                                + ' Invalid species id for initial condition. '
                                + 'It should be a valid index in the list of species.' 
                            )
                        f.write(f"{ins.speciesid+1}\n")

                    # And the concentration distribution
                    f.write(ins.concentration.get_file_entry())

        # Done
        return
