'''
Configuration of MODPATH-RW initial conditions
'''

# python
import numpy as np
from enum import Enum

# flopy
from flopy.pakbase import Package
from flopy.utils import Util3d

# local 
from .utils import multipackage


class icKindOption(Enum): 
    '''
    Enumerate formats for timeseries output option
    '''
    globalmass = 0
    localmass  = 1


class particlesDistributionOption(Enum): 
    '''
    Enumerate formats for timeseries output option
    '''
    equ         = 0 
    equispaced  = 0
    qra         = 1 
    quasirandom = 1
    ran         = 2
    random      = 2


class ModpathRWIc( Package ):
    """
    MODPATH-RW Initial Condition Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class: ModpathRW) to which
        this package will be added.
    kind  : int or str 
        The kind of initial condition specification. All formats load concentration 
        with u3d. Allowed values are: 
          * 0 ('globalmass'): concentration is transformed into a distribution
            of particles, by considering the cell volume, porosity and retardation factor.
            Distribution is obtained by transforming the total mass inside the
            cell into a number of particles estimated by means of particlesmass,
            distributed according to the cell aspect ratio, while using equispaced 
            and quasi-random particle distributions. Particles mass is unique for all 
            particles, recomputed in order to impose consistency between the total mass
            contained in the distribution and the estimated with particles.
          * 1 ('localmass'): similar to the previous alternative, particlesmass is employed 
            to estimate a number of particles per cell. The difference in this case, 
            is that once a number of particles is determined, particles mass is recomputed
            by establishing mass consistency at a cell level. This means that particles 
            from different cells may have different mass. Global mass consistency is 
            automatically satisfied.
    particlesmass  : float
        The scale of particles mass used for transforming the concentration data into 
        a distribution of particles. Determines resolution. 
    concentration : list, ndarray, int or float
        Interpreted as a resident initial concentration. It will be written with Util3d, 
        meaning that if the input is multidimensional, it should be consistent with 
        flow model dimensions. 
    particlesdist : int or str 
        The protocol for placing the particles inside each cell. Allowed values are: 
          * 0 ('equ','equispaced'): particles are placed following a uniform template
            determined from the estimated number of particles and the cell aspect 
            ratio. Particles are equispaced, and at half distance from the cell walls. 
          * 1 ('qra','quasirandom'): the equispaced template is perturbated based on
            a random number following a uniform distribution U[-1,1], proportional 
            to half the particle spacing.
          * 2 ('ran','random'): particles are placed randomly inside the cell. 
    drape : int 
        Determines the treatment of particles placed on initially dry cells. 
        Allowed values are: 
          * 0 : particles are unreleased. 
          * 1 : particles are vertically transfered to the next uppermost active cell. 
        notes:
          The drape option only applies for MODFLOW models solved with the standard formulation. 
          Models solved with Newton-Raphson are not affected by this parameter.
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
        particlesdist = 0   ,
        drape         = 0   ,
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

        # Initial condition kind
        if ( not isinstance( kind, (int,str) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for initial condition kind."
                f" It can be specified as int or str, but {str(type(kind))}"
                f" was given." 
            )
        if ( isinstance( kind, int ) ):
            # Determines format for establishing mass consistency
            if (kind not in [0,1]): 
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid initial condition kind parameter {str(kind)}."
                    f" Allowed values are 0 (globalmass) or 1 (localmass)."
                )
            self.kind = kind
        elif ( isinstance( kind, str ) ):
            try:
                self.kind = icKindOption[kind.lower()].value
            except:
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid initial condition kind {str(kind)}."
                    f" The allowed values are globalmass or localmass." 
                )

        # The particles mass
        if ( not isinstance( particlesmass, (int,float) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for particlesmass. It should be int/float,"
                f" but {str(type(particlesmass))} was given."
            )
        if (particlesmass <= 0):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid particle mass for initial condition."
                f" It should be a positive value." 
            )
        self.particlesmass = particlesmass

        # Particles distribution format 
        if ( not isinstance( particlesdist, (int,str) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for particlesdist parameter."
                f" It can be specified as int or str, but {str(type(particlesdist))}"
                f" was given." 
            )
        if ( isinstance( particlesdist, int ) ):
            # Determines format for establishing mass consistency
            if (particlesdist not in [0,1,2]): 
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid particlesdist parameter {str(particlesdist)}."
                    f" Allowed values are 0 (equispaced), 1 (quasirandom) "
                    f" or 2 (random)."
                )
            self.particlesdist = particlesdist
        elif ( isinstance( particlesdist, str ) ):
            try:
                self.particlesdist = particlesDistributionOption[particlesdist.lower()].value
            except:
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid particlesdist parameter {str(particlesdist)}."
                    f" The allowed values are equ/equispaced, qra/quasirandom"
                    f" or ran/random."
                )

        # Drape option 
        if ( drape is not None ):
            if not isinstance( drape, int ): 
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for drape. It should be integer,"
                    f" but {str(type(drape))} was given."
                )
            if (drape not in [0,1]):
                raise ValueError(
                    f"{self.__class__.__name__}:" 
                    f" Invalid value for drape. It should be"
                    f" 0 or 1."
                )
        else:
            drape = 0
        self.drape = drape

        # A species id, will only be written 
        # if particlesmassoption == 2
        if ( speciesid is not None ):
            if not isinstance( speciesid, int ): 
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for speciesid. It should be integer,"
                    f" but {str(type(speciesid))} was given."
                )
            if (speciesid < 0):
                raise ValueError(
                    f"{self.__class__.__name__}:" 
                    f" Invalid value for speciesid. It should be"
                    f" a valid zero-based index of species."
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
        ):
            # TypeError if any other unsupported type
            if ( not isinstance( concentration, (list,float,int,np.ndarray) ) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type {str(type(concentration))}"
                    f" for concentration. It should be list, np.ndarray or a single numeric value."
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
                        f"{self.__class__.__name__}:"
                        f" Invalid dtype {str(concentration.dtype)}" 
                        f" for concentration. It should be int or float."
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
                    f"{self.__class__.__name__}:"
                    f" Error while initializing distributed variable concentration."
                    f" Is the input shape consistent with flow model dimensions ?"
                )

        elif ( (self.kind in [0,1]) and (concentration is None) ):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid initial concentration. None was given."
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

                # Kind/format of initial condition, particles distribution and drape
                f.write(f"{ins.kind}   {ins.particlesdist}   {ins.drape}\n")

                # 0,1: resident concentration array  
                if ( ins.kind in [0,1] ):

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
                                f"{self.__class__.__name__}:"
                                f" Invalid species id for initial condition."
                                f" It should be a valid index in the list of species."
                            )
                        f.write(f"{ins.speciesid+1}\n")

                    # And the concentration distribution
                    f.write(ins.concentration.get_file_entry())

        # Done
        return
