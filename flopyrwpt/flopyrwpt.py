'''
Implement classes and methods for configuring mpathrwpt sims with flopy
'''


import flopy 
from flopy.pakbase import Package
from flopy.utils import Util2d, Util3d
from enum import Enum
import numpy as np
import pandas as pd
import os
from flopy.modpath.mp7particlegroup import (
    ParticleGroup,
    ParticleGroupLRCTemplate,
    ParticleGroupNodeTemplate,
)


from flopy.modpath.mp7particledata import (
    ParticleData,
)



class ModpathRW( flopy.modpath.Modpath7 ):
    '''
    MODPATH-RW Main class 
    Simply extends flopy.modpath.Modpath7
    '''
    def __init__(self, *args, **kwargs):
        # Call parent constructor
        super().__init__(*args,**kwargs)


class ModpathRWBas( flopy.modpath.Modpath7Bas ):
    '''
    MODPATH-RW Main class 
    Simply extends flopy.modpath.Modpath7
    '''
    def __init__(self, *args, **kwargs):
        # Call parent constructor
        super().__init__(*args,**kwargs)



class ModpathRWDsp( Package ):
    """
    MODPATH RW Dispersion Package Class.
    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    longitudinal : float or array of floats (nlay, nrow, ncol)
        The longitudinal dispersivity array (the default is XX).
    transverse : float or array of floats (nlay, nrow, ncol)
        The transverse dispersivity array (the default is XX).
    molecular : float 
        Molecular diffusion corrected by tortuosity effects
        The transverse dispersivity array (the default is 0).
    advection : str
        Advection model ot be used for advective component in RWPT integration.
        Could be set to 'eulerian' or 'exponential' (the default is 'eulerian').
    dimensions : int
        Define whether RWPT displacement should be done in 2D or 3D. 
        For 2D models, layers should be related to z direction and specify 2.
        (the default is 3).
    timestep : str
        Define method for computing particles time step. Can be courant, peclet or min. 
        The latter selects the minimum estimated between courant and peclet. For each 
        kind, it is expected that courant and/or peclet parameters are specified.
        (the default is courant with courant = 0.1)
    courant  : float
        Courant number that will be used at each cell for computing time step 
    peclet   :  float 
        Peclet number occupied for computing time step
    extension : str, optional
        File extension (default is 'dispersion').
    """

    def __init__(
        self,
        model,
        displong          = 0.0  , # xx
        disptranh         = 0.0  , # yy
        disptranv         = 0.0  , # zz
        diffaqueous       = 0.0  ,
        diffeff           = 0.0  , # corrected by tortuosity
        modelkind         = None , # linear, nonlinear
        betatransverse    = 0.5  , 
        betalongitudinal  = 1    ,
        mediumdelta       = 5    , 
        mediumdsize       = 1    ,
        extension         = 'dsp',
    ):

        pass



class ModpathRWOptions( Package ):
    """
    MODPATH-RW Options Package Class.
    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    longitudinal : float or array of floats (nlay, nrow, ncol)
        The longitudinal dispersivity array (the default is XX).
    transverse : float or array of floats (nlay, nrow, ncol)
        The transverse dispersivity array (the default is XX).
    molecular : float 
        Molecular diffusion corrected by tortuosity effects
        The transverse dispersivity array (the default is 0).
    advection : str
        Advection model ot be used for advective component in RWPT integration.
        Could be set to 'eulerian' or 'exponential' (the default is 'eulerian').
    dimensions : int
        Define whether RWPT displacement should be done in 2D or 3D. 
        For 2D models, layers should be related to z direction and specify 2.
        (the default is 3).
    timestep : str
        Define method for computing particles time step. Can be courant, peclet or min. 
        The latter selects the minimum estimated between courant and peclet. For each 
        kind, it is expected that courant and/or peclet parameters are specified.
        (the default is courant with courant = 0.1)
    courant  : float
        Courant number that will be used at each cell for computing time step 
    peclet   :  float 
        Peclet number occupied for computing time step
    extension : str, optional
        File extension (default is 'dispersion').
    """

    def __init__(
        self,
        model,
        dimensionmask     = [1,1,1],
        timestep          = 'courant',
        courant           = 0.1,
        ctdisp            = 0.1,
        advection         = 'eulerian',
    ):
    

        unitnumber = model.next_unit()
        super().__init__(model, extension, "RWOPTS", unitnumber)

        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])
        self.model = model
        self.shape3d = shape3d

        
        self.dimensionmask = dimensionmask




        # Minor checkings for given parameters 
        if ( advection not in ['eulerian','exponential'] ):
            raise Exception('flopyrwpt.py: advection model ' + advection + ' is not valid. Use eulerian or exponential')
        self.advection  = advection

        if ( dimensions not in [2,3] ):
            raise Exception('flopyrwpt.py: dimensions ' + str(dimensions) + ' is not valid. Use 2 or 3')
        self.dimensions = dimensions

        if ( timestep not in ['courant', 'peclet', 'min_adv_disp'] ):
            raise Exception('flopyrwpt.py: time step selection model ' + timestep + ' is not valid. courant, peclet, min_adv_disp')
        self.timestep  = timestep


        # Trust
        self.courant = courant
        self.peclet  = peclet 
        self.ctdisp  = ctdisp 


        self.parent.add_package(self)


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


        # Write initial conditions (IC)
        if self.nics > 0:
            # How many ics
            f.write(f"{self.nics}\n")

            # Loop over ic, notice pg
            for idic, ic in enumerate(self.initialconditions):
                # Write the initial condition
                ic.write(
                  f=f,
                  particlesmassoption=self.particlesmassoption,
                  solutesoption=self.solutesoption,
                )
        else:
            f.write(f"0\n")


        # Write prescribed flux boundaries
        if self.npfs > 0:
            # How many ics
            f.write(f"{self.npfs}\n")

            # Loop over pf
            for idpf, pf in enumerate(self.fluxconditions):
                # Write the flux condition 
                pf.write(
                  f=f,
                  particlesmassoption=self.particlesmassoption,
                  solutesoption=self.solutesoption,
                )
        else:
            f.write(f"0\n")


        # Write icbound
        f.write(self.icbound.get_file_entry())


        # Depending on the solutes option 
        # is how to continue writing the file
        # self.solutesoption is assigned in ModpathRwpSim
        if self.solutesoption == 0:
            # Single virtual solute, all 
            # pgroups are displaced with the 
            # same dispersion properties

            # Select model dispersivity
            if self.modelkind == 1:

                # Write modelkind
                f.write(f"{self.modelkind}\n")

                # Write molecular diffusion
                f.write(f"{self.molecular:.16f}\n")

                # Write dispersivities
                f.write(self.longitudinal.get_file_entry())
                f.write(self.transverse.get_file_entry())

            elif self.modelkind == 2:

                # Write modelkind
                f.write(f"{self.modelkind}\n")

                # Write nonlinear beta exponents
                f.write(f"{self.betalongitudinal:.10f}\n") 
                f.write(f"{self.betatransverse:.10f}\n") 
                f.write(f"{self.mediumdelta:.10f}\n")
                #f.write(self.betalongitudinal.get_file_entry())
                #f.write(self.betatransverse.get_file_entry())
                #f.write(self.mediumdelta.get_file_entry())
                f.write(self.mediumdsize.get_file_entry())

                # Write aqueous molecular diffusion
                f.write(f"{self.molecular:.16f}\n")


        elif self.solutesoption ==1:

            f.write(f"{len(self.solutes)}\n")

            # Multiple solutes
            for ns, sol in enumerate(self.solutes):
                sol.write(f=f, particlesmassoption=self.particlesmassoption)
        

        # Write time step selection method
        if self.timestep is not None:
            if self.timestep == 'peclet': 
                f.write(f"CONSTANT_PE\n")
                f.write(f"{self.peclet:.10f}\n")  
            if self.timestep == 'courant': 
                f.write(f"CONSTANT_CU\n")  
                f.write(f"{self.courant:.10f}\n")  
            if self.timestep == 'min_adv_disp':
                # Review this format, it may require both
                f.write(f"MIN_ADV_DISP\n")  
                f.write(f"{self.courant:.10f}\n")  
                f.write(f"{self.ctdisp:.10f}\n")  
                #f.write(f"{self.courant:20d}") #? 

        # Write advection model
        f.write(f"{self.advection.upper():20s}\n")  

        # Write number of dimensions for RWPT displacement
        f.write(f"{self.dimensions:1d}D\n")  


        # And close
        f.close()


class ModpathRwptDispersion( Package ):
    """
    MODPATH RWPT Dispersion Package Class.
    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    longitudinal : float or array of floats (nlay, nrow, ncol)
        The longitudinal dispersivity array (the default is XX).
    transverse : float or array of floats (nlay, nrow, ncol)
        The transverse dispersivity array (the default is XX).
    molecular : float 
        Molecular diffusion corrected by tortuosity effects
        The transverse dispersivity array (the default is 0).
    advection : str
        Advection model ot be used for advective component in RWPT integration.
        Could be set to 'eulerian' or 'exponential' (the default is 'eulerian').
    dimensions : int
        Define whether RWPT displacement should be done in 2D or 3D. 
        For 2D models, layers should be related to z direction and specify 2.
        (the default is 3).
    timestep : str
        Define method for computing particles time step. Can be courant, peclet or min. 
        The latter selects the minimum estimated between courant and peclet. For each 
        kind, it is expected that courant and/or peclet parameters are specified.
        (the default is courant with courant = 0.1)
    courant  : float
        Courant number that will be used at each cell for computing time step 
    peclet   :  float 
        Peclet number occupied for computing time step
    extension : str, optional
        File extension (default is 'dispersion').
    """

    def __init__(
        self,
        model,
        longitudinal      = 1,
        transverse        = 0.1,
        molecular         = 0.0,
        advection         = 'eulerian',
        dimensions        = 3,
        timestep          = 'courant',
        courant           = 0.1,
        peclet            = 100.0,
        ctdisp            = 0.1,
        extension         = 'dispersion',
        icbound           = None,
        modelkind         = None, # linear, nonlinear
        betatransverse    = 0.5, 
        betalongitudinal  = 1,
        mediumdelta       = 5, 
        mediumdsize       = 1,
        solutes             = None, 
    ):
    

        unitnumber = model.next_unit()
        super().__init__(model, extension, "DISPERSION", unitnumber)

        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])
        self.model = model
        self.shape3d = shape3d

        # Needs some health checking
        self.solutes = solutes 

        # Select modelkind
        if modelkind is not None: 
            if modelkind not in ('linear', 'nonlinear'):
                raise Exception( 'flopyrwpt.py: modelkind ' + modelkind + ' is not allowed: linear or nonlinear' )
            if modelkind == 'linear':
                self.modelkind = 1
            elif modelkind == 'nonlinear':
                self.modelkind = 2
        else:
            self.modelkind = 1

        # Assign NONLINEAR dispersion parameters
        self.betalongitudinal = betalongitudinal
        self.betatransverse = betatransverse
        self.mediumdelta = mediumdelta
        ## Assign NONLINEAR dispersion parameters
        #self.betalongitudinal = Util3d(
        #    model,
        #    shape3d,
        #    np.float32,
        #    betalongitudinal,
        #    name="BETALONGITUDINAL",
        #    locat=self.unit_number[0],
        #)
        #self.betatransverse= Util3d(
        #    model,
        #    shape3d,
        #    np.float32,
        #    betatransverse,
        #    name="BETATRANSVERSE",
        #    locat=self.unit_number[0],
        #)
        #self.mediumdelta = Util3d(
        #    model,
        #    shape3d,
        #    np.float32,
        #    mediumdelta,
        #    name="MEDIUMDELTA",
        #    locat=self.unit_number[0],
        #)
        self.mediumdsize = Util3d(
            model,
            shape3d,
            np.float32,
            mediumdsize,
            name="MEDIUMDSIZE",
            locat=self.unit_number[0],
        )


        # Assign LINEAR dispersion parameters
        self.longitudinal = Util3d(
            model,
            shape3d,
            np.float32,
            longitudinal,
            name="LONGITUDINAL",
            locat=self.unit_number[0],
        )
        self.transverse= Util3d(
            model,
            shape3d,
            np.float32,
            transverse,
            name="TRANSVERSE",
            locat=self.unit_number[0],
        )
        self.molecular  = molecular


        # ICBOUND
        if icbound is not None:
            self.icbound= Util3d(
                model,
                shape3d,
                np.int32,
                icbound,
                name="ICBOUND",
                locat=self.unit_number[0],
            )
        else:
            self.icbound= Util3d(
                model,
                shape3d,
                np.int32,
                0,
                name="ICBOUND",
                locat=self.unit_number[0],
            )


        # Minor checkings for given parameters 
        if ( advection not in ['eulerian','exponential'] ):
            raise Exception('flopyrwpt.py: advection model ' + advection + ' is not valid. Use eulerian or exponential')
        self.advection  = advection

        if ( dimensions not in [2,3] ):
            raise Exception('flopyrwpt.py: dimensions ' + str(dimensions) + ' is not valid. Use 2 or 3')
        self.dimensions = dimensions

        if ( timestep not in ['courant', 'peclet', 'min_adv_disp'] ):
            raise Exception('flopyrwpt.py: time step selection model ' + timestep + ' is not valid. courant, peclet, min_adv_disp')
        self.timestep  = timestep


        # Trust
        self.courant = courant
        self.peclet  = peclet 
        self.ctdisp  = ctdisp 


        self.parent.add_package(self)


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


        # Write initial conditions (IC)
        if self.nics > 0:
            # How many ics
            f.write(f"{self.nics}\n")

            # Loop over ic, notice pg
            for idic, ic in enumerate(self.initialconditions):
                # Write the initial condition
                ic.write(
                  f=f,
                  particlesmassoption=self.particlesmassoption,
                  solutesoption=self.solutesoption,
                )
        else:
            f.write(f"0\n")


        # Write prescribed flux boundaries
        if self.npfs > 0:
            # How many ics
            f.write(f"{self.npfs}\n")

            # Loop over pf
            for idpf, pf in enumerate(self.fluxconditions):
                # Write the flux condition 
                pf.write(
                  f=f,
                  particlesmassoption=self.particlesmassoption,
                  solutesoption=self.solutesoption,
                )
        else:
            f.write(f"0\n")


        # Write icbound
        f.write(self.icbound.get_file_entry())


        # Depending on the solutes option 
        # is how to continue writing the file
        # self.solutesoption is assigned in ModpathRwpSim
        if self.solutesoption == 0:
            # Single virtual solute, all 
            # pgroups are displaced with the 
            # same dispersion properties

            # Select model dispersivity
            if self.modelkind == 1:

                # Write modelkind
                f.write(f"{self.modelkind}\n")

                # Write molecular diffusion
                f.write(f"{self.molecular:.16f}\n")

                # Write dispersivities
                f.write(self.longitudinal.get_file_entry())
                f.write(self.transverse.get_file_entry())

            elif self.modelkind == 2:

                # Write modelkind
                f.write(f"{self.modelkind}\n")

                # Write nonlinear beta exponents
                f.write(f"{self.betalongitudinal:.10f}\n") 
                f.write(f"{self.betatransverse:.10f}\n") 
                f.write(f"{self.mediumdelta:.10f}\n")
                #f.write(self.betalongitudinal.get_file_entry())
                #f.write(self.betatransverse.get_file_entry())
                #f.write(self.mediumdelta.get_file_entry())
                f.write(self.mediumdsize.get_file_entry())

                # Write aqueous molecular diffusion
                f.write(f"{self.molecular:.16f}\n")


        elif self.solutesoption ==1:

            f.write(f"{len(self.solutes)}\n")

            # Multiple solutes
            for ns, sol in enumerate(self.solutes):
                sol.write(f=f, particlesmassoption=self.particlesmassoption)
        

        # Write time step selection method
        if self.timestep is not None:
            if self.timestep == 'peclet': 
                f.write(f"CONSTANT_PE\n")
                f.write(f"{self.peclet:.10f}\n")  
            if self.timestep == 'courant': 
                f.write(f"CONSTANT_CU\n")  
                f.write(f"{self.courant:.10f}\n")  
            if self.timestep == 'min_adv_disp':
                # Review this format, it may require both
                f.write(f"MIN_ADV_DISP\n")  
                f.write(f"{self.courant:.10f}\n")  
                f.write(f"{self.ctdisp:.10f}\n")  
                #f.write(f"{self.courant:20d}") #? 

        # Write advection model
        f.write(f"{self.advection.upper():20s}\n")  

        # Write number of dimensions for RWPT displacement
        f.write(f"{self.dimensions:1d}D\n")  


        # And close
        f.close()




class ModpathRwptReconstruction( Package ):
    """
    MODPATH RWPT Dispersion Package Class.
    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    binsize : list [sx,sy,sz]
        List with cell size employed for reconstruction process 
    domainsize : list [sx,sy,sz]
        List with domain sizes. These are used for computing reconstruction 
        grid dimensions using given binsize.
    domainorigin : list [ox,oy,oz]
        List with domain origin. 
    noptloops: int
        The number of optimization loops for determinin smoothing conditions
    outputfilename: str
        Filename that will be used for writing reconstruction output. Contains 
        cell id, density and both time step and particle group identification.
    extension : str, optional
        File extension (default is 'dispersion').
    """

    def __init__(
        self,
        model,
        binsize              = [1,1,1],
        domainsize           = [1,1,1],
        domainorigin         = [0,0,0],
        kerneldatabase       = False, 
        kerneldatabaseparams = [1,0.1,30], # minh/lambda,deltahlambda,maxhlambda
        noptloops            = 10, 
        outputfilename       = 'gpkde.output',
        extension            = 'gpkde',
        skiptimeserieswriter = False,
        minhlambda           = None,
        maxhlambda           = None,
        deltahlambda         = None,
    ):

        unitnumber = model.next_unit()

        super().__init__(model, extension, "GPKDE", unitnumber)

        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])

        self.binsize      = np.array(binsize).astype(np.float32)
        self.domainsize   = np.array(domainsize).astype(np.float32)
        self.domainorigin = np.array(domainorigin).astype(np.float32)
        self.noptloops    = noptloops

        self.kerneldatabase       = kerneldatabase       
        self.kerneldatabaseparams = kerneldatabaseparams 

        self.skiptimeserieswriter = skiptimeserieswriter

        if outputfilename is not None:
            self.outputfilename = outputfilename
        else:
            raise Exception('requires outputfilename for reconstruction output')

        self.parent.add_package(self)

        

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

        # Output filename
        f.write(f"{self.outputfilename}\n")

        # Skip timeseries writer 
        if self.skiptimeserieswriter:
            f.write(f"1\n") # 1 true, skip TimeseriesWriter
        else:
            f.write(f"0\n") # 0 false, write timeseries

        # Domain origin
        for idb, b in enumerate(self.domainorigin):
            if idb == len(self.domainorigin)-1: 
                f.write(f"{b:10f}\n")
            else:
                f.write(f"{b:10f}\t")

        # Domain sizes
        for idb, b in enumerate(self.domainsize):
            if idb == len(self.domainsize)-1: 
                f.write(f"{b:10f}\n")
            else:
                f.write(f"{b:10f}\t")

        # Bin sizes
        for idb, b in enumerate(self.binsize):
            if idb == len(self.binsize)-1: 
                f.write(f"{b:10f}\n")
            else:
                f.write(f"{b:10f}\t")

        # Number of optimization loops 
        f.write(f"{self.noptloops:10d}\n")

        # Kernel database reconstruction or direct without database
        # 1: kerneldatabase
        # 0: bruteforcedensity
        if self.kerneldatabase:
            # Kernel database reconstruction
            f.write(f"1\n") # 1 for id into fortran

            # Database params: minh/lambda, daltah/lambda, maxh/lambda
            for idb, b in enumerate(self.kerneldatabaseparams):
                if idb == len(self.binsize)-1: 
                    f.write(f"{b:10f}\n")
                else:
                    f.write(f"{b:10f\t}")
        else:
            # Brute force reconstruction
            f.write(f"0\n") # 0 for id into fortran
            
            # kernel params: minh/lambda, maxh/lambda
            for idb, b in enumerate(self.kerneldatabaseparams):
                if idb == 0:
                    f.write(f"{b:16f}\t")
                if idb == 2: 
                    f.write(f"{b:16f}\n")


        # And close
        f.close()





def count_instances(method):
    def wrapper(self, *args, **kw):
        self.__class__.COUNTER += 1
        return method(self, *args, **kw)
    return wrapper


class ModpathRwptObs( Package ):
    """
    MODPATH RWPT Observation Package Class.
    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    """


    COUNTER = 0

    @count_instances
    def __init__(
        self,
        model,
        id             = None,
        kind           = 1,
        celloption     = 1,
        cells          = None,
        timeoption     = 1,
        structured     = True, 
        basefilename   = 'mpathrwobs_',
        filename       = None,
        extension      = '.obs',
    ):
        
        unitnumber = model.next_unit()
        super().__init__(model, extension, "OBSCELLS", unitnumber)

        self.kind = kind
        self.timeoption = timeoption
        self.structured = structured

        if ( celloption not in [1,2] ):
            raise ValueError('Invalid celloption ',
                    celloption, '. Allowed values are 1 (list of cellids) or 2 (modelgrid like array)')
        self.celloption = celloption

        # Write a list of cellids
        if self.celloption == 1:
            if cells is None:
                self.cells = []
            else:
                if not isinstance( cells, (list,np.ndarray,np.array) ):
                    raise Exception( 'Cells parameter should be a list or numpy array')
                # Maybe some sanity check about data structure or the same 
                # used for partlocs
                self.cells = cells

        # Write obs cells as 3D array
        if self.celloption == 2: 

            # This was already done right?
            shape = model.shape
            if len(shape) == 3:
                shape3d = shape
            elif len(shape) == 2:
                shape3d = (shape[0], 1, shape[1])
            else:
                shape3d = (1, 1, shape[0])
            self.model   = model
            self.shape3d = shape3d

            if cells is None:
                self.cells = []
            else:
                if not isinstance( cells, (list,np.ndarray) ):
                    raise Exception( 'Cells parameter should be a list or numpy array')
                self.cells = Util3d(
                    model,
                    shape3d,
                    np.int32,
                    cells,
                    name="OBSCELLS",
                    locat=self.unit_number[0],
                )

        # Define obs id
        # It could be more useful some kind of random 
        # integer
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER


        self.filename = basefilename+str(self.id)+extension


        return


    def write(self, f=None):
        """
        Write the package file
        Parameters
        ----------
        Returns
        -------
        None
        """


        if f is None:
            raise Exception('ModpathRwptObs: requires file pointer f. Is None.')

        # validate that a valid file object was passed
        if not hasattr(f, "write"):
            raise ValueError(
                "{}: cannot write data for template without passing a valid "
                "file object ({}) open for writing".format(self.name, f)
            )


        # Write obs id
        f.write(f"{self.id}\n")

        # Write obs filename
        f.write(f"{self.filename}\n")

        # Write the obs kind
        f.write(f"{self.kind}\n")

        # Write the celloption param
        f.write(f"{self.celloption}\n")

        if self.celloption == 1:
            # Should write a cell number
            if len( self.cells ) == 0: 
                raise Exception('Observation cells is empty. Specify a list of cells for the observation id ', self.id)
            else:
                f.write(f"{len(self.cells)}\n")
                fmts = []
                if self.structured:
                    f.write(f"1\n") # To indicate structured
                    fmts.append("{:9d}") # lay
                    fmts.append("{:9d}") # row
                    fmts.append("{:9d}") # col
                else:
                    f.write(f"2\n") # To indicate cell ids
                    fmts.append("{:9d}") # cellid
                fmt = " " + " ".join(fmts) + "\n"
                for oc in self.cells:
                    woc = np.array(oc).astype(np.int32)+1 # Correct the zero-based indexes
                    f.write(fmt.format(*woc))

        elif self.celloption == 2:
            # Should write an array
            # with the distribution of the observation
            f.write(self.cells.get_file_entry())


        # Write timeoption params
        f.write(f"{self.timeoption}\n")
        if self.timeoption == 1:
            pass
        elif self.timeoption == 2:
            pass




        return


class ModpathRwptSolute( Package ):
    """
    MODPATH RWPT Solute Package Class.
    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    """

    COUNTER = 0

    @count_instances
    def __init__(
        self,
        model,
        id       = None,
        stringid = None,
        kind     = 1,
        dispmodel= 1,
        daqueous = 0.0,
        displong = None,
        disptranh= None,
        disptranv= None,
        pgroups  = None,
        extension='.sol',
    ):
        
        unitnumber = model.next_unit()
        super().__init__(model, extension, "SOLUTE", unitnumber)

        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])

        self.model = model
        self.shape3d = shape3d

        self.kind     = kind
        self.daqueous = daqueous
        
        self.dispmodel = dispmodel
       
        # Assign LINEAR dispersion parameters
        self.displong = Util3d(
            model,
            shape3d,
            np.float32,
            displong,
            name="LONGITUDINAL",
            locat=self.unit_number[0],
        )
        if(
           ( disptranh is None ) and 
           ( disptranv is not None ) ):
            self.disptranh = disptranv
        else:
            self.disptranh= Util3d(
                model,
                shape3d,
                np.float32,
                disptranh,
                name="TRANSVERSEH",
                locat=self.unit_number[0],
            )
        if(
           ( disptranv is None ) and 
           ( disptranh is not None ) ):
            self.disptranv = self.disptranh
        else:
            self.disptranv= Util3d(
                model,
                shape3d,
                np.float32,
                disptranv,
                name="TRANSVERSEV",
                locat=self.unit_number[0],
            )


        if not isinstance( pgroups, (list, np.ndarray)):
            raise TypeError('pgroups should be a list or array, given ', type(pgroups))
        self.pgroups = pgroups



        # Define obs id
        # It could be more useful some kind of random 
        # integer
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER

        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = 'S'+str(self.__class__.COUNTER)


        return


    def write(self, f=None, particlesmassoption=0):

        """
        Write the package file
        Parameters
        ----------
        Returns
        -------
        None
        """


        if f is None:
            raise Exception('ModpathRwptObs: requires file pointer f. Is None.')

        # validate that a valid file object was passed
        if not hasattr(f, "write"):
            raise ValueError(
                "{}: cannot write data for template without passing a valid "
                "file object ({}) open for writing".format(self.name, f)
            )


        # Write id
        f.write(f"{self.id}\n")

        # Write stringid
        f.write(f"{self.stringid}\n")

        # Need to write the related groups 
        # only if particlesmassoption not equal 2 
        if ( particlesmassoption != 2 ):
            # Write the number of groups
            f.write(f"{len(self.pgroups)}\n")

            for idpg, pg in enumerate(self.pgroups):
                # Write the pgroup id in the list of pgroups
                f.write(f"{pg+1}\n")

        if self.dispmodel == 1:

            # Write modelkind
            f.write(f"{self.dispmodel}\n")

            # Write effetive molecular diffusion
            f.write(f"{self.daqueous:.16f}\n")

            # Write dispersivities
            f.write(self.displong.get_file_entry())
            f.write(self.disptranh.get_file_entry())

        elif self.dispmodel == 2:
            pass


        return


class ModpathRwptIc( Package ):
    """
    MODPATH RWPT Initial Condition Package Class.
    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    concentration : ndarray
        Interpreted as resident concentration, that is, the product beween porosity and
        dissolved concentration. 
    """

    COUNTER = 0

    @count_instances
    def __init__(
        self,
        model,
        id            = None,
        stringid      = None,
        kind          = 1,
        mass          = 1,
        massformat    = 1,
        concentration = None,
        soluteid      = 1,
        extension     = '.ic',
    ):
        
        unitnumber = model.next_unit()
        super().__init__(model, extension, "RWIC", unitnumber)

        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])

        self.model = model
        self.shape3d = shape3d

        # Determines how to write the initial condition
        if (kind not in [1]): 
            raise ValueError('ModpathRWIc: Invalid value for kind, should be 1. Given ', str(kind) )
        self.kind = kind

        self.mass = mass
        self.massformat = massformat

        self.soluteid = soluteid


        # Generic 
        # id
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER
        # stringid 
        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = 'IC'+str(self.__class__.COUNTER)

        if (
            ( concentration is not None ) and 
            ( self.kind == 1 )
        ):
            if not isinstance( concentration, np.ndarray ):
                raise TypeError('ModpathRWIc: concentration should be an np.ndarray')
            self.concentration = Util3d(
                self.model,
                self.shape3d,
                np.float32,
                concentration,
                name=self.stringid,
                locat=self.unit_number[0],
            )
        elif ( (self.kind ==1) and (concentration is None) ):
            raise ValueError('ModpathRWIc: kind 1 requires concentration array, None given.')
        
        return


    def write(self, f=None, particlesmassoption=0, solutesoption=0):

        """
        Write the package file
        Parameters
        ----------
        Returns
        -------
        None
        """


        # validate that a valid file object was passed
        if not hasattr(f, "write"):
            raise ValueError(
                "{}: cannot write data for template without passing a valid "
                "file object ({}) open for writing".format(self.name, f)
            )

        # Write string id 
        f.write(f"{self.stringid}\n")

        # Kind/format of initial condition 
        f.write(f"{self.kind}\n")

        # 1: resident concentration array  
        if self.kind == 1:

            # Give particles mass 
            f.write(f"{self.mass:.10f}\n")

            # If solutes are being specified, do it
            if( (particlesmassoption == 2) or (solutesoption == 1)):
                f.write(f"{self.soluteid}\n")

            # And concentration distribution
            f.write(self.concentration.get_file_entry())

        return


class ModpathRwptFlux( Package ):
    """
    MODPATH RWPT Flux Package Class.
    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    """

    COUNTER = 0

    @count_instances
    def __init__(
        self,
        model,
        id            = None,
        stringid      = None,
        kind          = 1,
        mass          = 1,
        massformat    = 1,   # could be used eventually
        starttime     = None, 
        endtime       = None,
        dtrelease     = None,
        times         = None,
        concentration = None,
        cellid        = None,
        structured    = True,
        nparticles    = None,
        rowsubdiv     = 1,
        colsubdiv     = 1,
        laysubdiv     = 1,
        soluteid      = 1,
        filename      = None,
        extension     = '.flux',
    ):
        
        unitnumber = model.next_unit()
        super().__init__(model, extension, "RWFLUX", unitnumber)

        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])

        self.model = model
        self.shape3d = shape3d

        # Determines how to write the initial condition
        if (kind not in [1,2]): 
            raise ValueError(self.__class__.__name__,': Invalid value for kind, should be 1,2. Given ', str(kind) )
        # Concentration timeseries
        if kind == 2:
            if ( times is None ):
                raise ValueError(self.__class__.__name__,': Format 2 for injection requires times, but None was given.' )
            if not isinstance( times, (list,np.ndarray) ):
                raise ValueError(self.__class__.__name__,': Format 2 requires times as list or np.ndarray. Given ', type(times) )
            if ( concentration is None ):
                raise ValueError(self.__class__.__name__,': Format 2 for injection requires concentrations, but None was given.' )
            if not isinstance( concentration, (list,np.ndarray) ):
                raise ValueError(self.__class__.__name__,': Format 2 requires concentrations as list or np.ndarray. Given ', type(times) )
            if len(concentration) != len(times):
                raise ValueError(self.__class__.__name__,': Format 2 requires concentrations to have the same size than times. ' )
        self.kind = kind

        self.soluteid = soluteid
        self.structured = structured


        if cellid is None: 
            raise ValueError(self.__class__.__name__,': Requires cellid. Given ', str(cellid) )
        if isinstance(cellid, (list,np.ndarray)):
            if ( self.structured and (len(cellid)!=3) ):
                raise ValueError(self.__class__.__name__,': Structured cellid should follow [lay,row,col] format. Given ', str(cellid) )
            elif ( (not self.structured) and (len(cellid)!=1) ): 
                raise ValueError(self.__class__.__name__,': Structured cellid should follow [lay,row,col] format. Given ', str(cellid) )
        elif isinstance(cellid, (int)):
            if ( self.structured ): 
                raise ValueError(self.__class__.__name__,': Structured cellid requires [lay,row,col] format. Given ', str(cellid) )
        else:
            raise TypeError(self.__class__.__name__,': Cellid should be a list, ndarray or integer depending on structured parameter.' )
        self.cellid = cellid

        # Assign mass format
        if (massformat not in [1]): 
            raise ValueError(self.__class__.__name__,': Invalid value for massformat, should be 1. Given ', str(massformat) )
        self.massformat = massformat
        if self.massformat == 1:
            self.mass = mass
        #elif self.massformat == 2: # requires nparticles
        #    # Compute mass based on dtrelease, concentration, flow-rate
        #    if flowrate is None:
        #        raise ValueError(self.__class__.__name__,': Massformat 2 requires flow rate.' )
        #    mass = flowrate*concentration*dtrelease/nparticles 
        #    self.mass = mass
        self.massformat = massformat


        # Generic 
        # id
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER
        # stringid 
        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = 'PF'+str(self.__class__.COUNTER)

        # If kind 1, verify times
        if(self.kind == 1):
            if ( starttime is None ):
                raise ValueError(self.__class__.__name__,': Format 1 for injection requires starttime but None was given.' )
            if ( dtrelease is None ):
                raise ValueError(self.__class__.__name__,': Format 1 for injection requires dtrelease but None was given.' )
            if ( endtime is None ):
                raise ValueError(self.__class__.__name__,': Format 1 for injection requires endtime but None was given.' )
            self.starttime = starttime
            self.dtrelease = dtrelease
            self.endtime   = endtime
            self.releasecount = int(len(np.arange(starttime,endtime+dtrelease,dtrelease)))
    

        # Define a filename for writing the injection series
        if( (self.kind == 2) and (filename is None)):
            self.filename = self.stringid + '_timeseries.csv'



        self.rowsubdiv = rowsubdiv
        self.colsubdiv = colsubdiv
        self.laysubdiv = laysubdiv


        return


    def write(self, f=None, particlesmassoption=0, solutesoption=0):

        """
        Write the package file
        Parameters
        ----------
        Returns
        -------
        None
        """


        # validate that a valid file object was passed
        if not hasattr(f, "write"):
            raise ValueError(
                "{}: cannot write data for template without passing a valid "
                "file object ({}) open for writing".format(self.name, f)
            )

        # Write string id 
        f.write(f"{self.stringid}\n")

        # Kind/format of flux condition 
        f.write(f"{self.kind}\n")

        if self.kind == 1:
            # Cell
            if self.structured:
                f.write(f"1\n") # To indicate structured
                fmts.append("{:9d}") # lay
                fmts.append("{:9d}") # row
                fmts.append("{:9d}") # col
            else:
                f.write(f"2\n") # To indicate cell ids
                fmts.append("{:9d}") # cellid
            fmt = " " + " ".join(fmts) + "\n"
            woc = np.array(self.cellid).astype(np.int32)+1 # Correct the zero-based indexes
            f.write(fmt.format(*woc))

            # If particles mass is given, assumes
            # estimation performed via some 
            # method involving concentration

            # Give particles mass 
            f.write(f"{self.mass:.10f}\n")
            
            # If solutes are being specified, do it
            if( (particlesmassoption == 2) or (solutesoption == 1)):
                f.write(f"{self.soluteid}\n")

            fmts = []
            fmts.append("{:9d}")   # Release time count
            fmts.append("{:.10f}") # Initial release
            fmts.append("{:.10f}") # Release interval 
            fmt = " " + " ".join(fmts) + "\n"
            woc = [self.releasecount,self.starttime,self.dtrelease]
            f.write(fmt.format(*woc))

            fmt = " {} {} {}\n"
            line = fmt.format(
                self.colsubdiv,
                self.rowsubdiv,
                self.laysubdiv,
            )
            f.write(line)
            

        # A injection series
        if self.kind == 2:
            # Write a cell number
            fmts = []
            if self.structured:
                f.write(f"1\n") # To indicate structured
                fmts.append("{:9d}") # lay
                fmts.append("{:9d}") # row
                fmts.append("{:9d}") # col
            else:
                f.write(f"2\n") # To indicate cell ids
                fmts.append("{:9d}") # cellid
            fmt = " " + " ".join(fmts) + "\n"
            woc = np.array(self.cellid).astype(np.int32)+1 # Correct the zero-based indexes
            f.write(fmt.format(*woc))

            # Give particles mass 
            f.write(f"{self.mass:.10f}\n")

            # If solutes are being specified, do it
            if( (particlesmassoption == 2) or (solutesoption == 1)):
                f.write(f"{self.soluteid}\n")

            # Write length of the timeseries
            f.write(f"{len(self.times)}\n")

            # Write the timeseries file
            df         = pd.DataFrame()
            df['time'] = self.times
            df['conc'] = self.concentration
            df.to_csv(os.path.join(
                self.model.model_ws, self.filename),
                header=None,
                index =False,
                sep   =' ', 
                float_format="%.10f"
            )
            # And its name
            f.write(f"{self.filename}\n")


        return



# Update particle groups classes
mp7ParticleGroup = ParticleGroup
mp7ParticleGroupLRCTemplate = ParticleGroupLRCTemplate
mp7ParticleGroupNodeTemplate = ParticleGroupNodeTemplate
class ParticleGroup( mp7ParticleGroup ): 

    def __init__( self, *args, mass=1.0, solute=0, **kwargs ):
        # Call parent constructor
        super().__init__(*args,**kwargs) 
        
        if mass != 0:
            self.mass = mass
        else:
            raise ValueError('Mass for the particle group is zero')

        self.solute = solute

    def write( self, fp=None, ws=".", mass=False, solute=False): 
        # Call base class write method to write common data
        super().write(fp, ws)

        if mass:
            # Write the particle mass
            fp.write(f"{self.mass:.16f}\n")

        if solute:
            # Write the solute id 
            fp.write(f"{self.solute:9df}\n")
        return

class ParticleGroupLRCTemplate( mp7ParticleGroupLRCTemplate ): 

    def __init__( self, *args, mass=1.0, solute=0, **kwargs ):
        # Call parent constructor
        super().__init__(*args,**kwargs) 
        
        if mass != 0:
            self.mass = mass
        else:
            raise ValueError('Mass for the particle group is zero')

        self.solute = solute

    def write( self, fp=None, ws=".", mass=False, solute=False): 
        # Call base class write method to write common data
        super().write(fp, ws)

        if mass:
            # Write the particle mass
            fp.write(f"{self.mass:.16f}\n")

        if solute:
            # Write the solute id 
            fp.write(f"{self.solute:9df}\n")
        return

class ParticleGroupNodeTemplate( mp7ParticleGroupNodeTemplate ): 

    def __init__( self, *args, mass=1.0, solute=0, **kwargs ):
        # Call parent constructor
        super().__init__(*args,**kwargs) 
        
        if mass != 0:
            self.mass = mass
        else:
            raise ValueError('Mass for the particle group is zero')

        self.solute = solute

    def write( self, fp=None, ws=".", mass=False, solute=False): 
        # Call base class write method to write common data
        super().write(fp, ws)

        if mass:
            # Write the particle mass
            fp.write(f"{self.mass:.16f}\n")
        if solute:
            # Write the solute id 
            fp.write(f"{self.solute:9d}\n")

        return


# Updates enumeration of available simulations
import flopy.modpath.mp7sim as mp7
class simType(Enum):
    """
    Enumeration of different simulation types
    """
    
    endpoint       = 1
    pathline       = 2
    timeseries     = 3
    combined       = 4
    rwpttimeseries = 5
    rwptcombined   = 6
    rwptendpoint   = 7
mp7.simType = simType



# Class for RWPT sims
class ModpathRwptSim( flopy.modpath.Modpath7Sim ):
    '''
    Should extend mp7sim to create mp7rwptsim
    '''

    def __init__( self, *args,
            timeseriesoutputoption=0,
            particlesmassoption=0,
            solutesoption=0,
            reconstruction=False,
            dispersionfilename=None,
            reconstructionfilename=None,
            observations=None,
            solutes=None,
            initialconditions=None,
            fluxconditions=None,
            **kwargs
        ):


        # Call parent constructor
        super().__init__(*args,**kwargs, extension='mprw') 


        # Override interpretation of particlegroups, allows None
        if 'particlegroups' not in kwargs.keys():
            particlegroups = None
        self.particlegroups = particlegroups


        # New options
        if (timeseriesoutputoption not in [0,1]):
            raise ValueError('Timeseries output option should be 0 or 1. Given :', str(timeseriesoutputoption) )
        self.timeseriesoutputoption = timeseriesoutputoption
        if (particlesmassoption not in [0,1,2]):
            raise ValueError('Particles mass option should be 0, 1 or 2. Given :', str(particlesmassoption) )
        self.particlesmassoption = particlesmassoption
        if (solutesoption not in [0,1]):
            raise ValueError('Solutes option should be 0 or 1. Given :', str(solutesoption) )
        self.solutesoption = solutesoption


        # Initial conditions
        if initialconditions is not None:
            if isinstance(initialconditions,list):
                for pic in initialconditions:
                    if not isinstance(pic, ModpathRwptIc):
                        raise TypeError('Object in list is type ', type(pic), '. Expected ModpathRwptIc.')
                # Survived so continue
                self.nics = len(initialconditions)
                self.initialconditions = initialconditions
            elif isinstance( initialconditions, ModpathRwptIc ):
                self.nics = 1
                self.initialconditions = [initialconditions]
            else:
                raise TypeError('Initial conditions argument should be of type list or ModpathRwptIc. ', type(initialconditions), ' given.')
        else:
            self.nics = 0
            self.initialconditions = None


        # Prescribed flux conditions
        if fluxconditions is not None:
            if isinstance(observations,list):
                for pic in fluxconditions:
                    if not isinstance(pic, ModpathRwptFlux):
                        raise TypeError('Object in list is type ', type(pic), '. Expected ModpathRwptFlux.')
                # Survived so continue
                self.npfs = len(fluxconditions)
                self.fluxconditions = fluxconditions
            elif isinstance( fluxconditions, ModpathRwptFlux ):
                self.npfs = 1
                self.fluxconditions = [fluxconditions]
            else:
                raise TypeError('Initial conditions argument should be of type list or ModpathRwptFlux. ', type(fluxconditions), ' given.')
        else:
            self.npfs = 0
            self.fluxconditions = None


        # Inform dispersion package about 
        # particlesmassoption
        # solutesoption
        if self.simulationtype > 4: 
            disp = self.parent.get_package('DISPERSION')
            if disp is None:
                raise Exception('Requires a dispersion package')
            disp.particlesmassoption = self.particlesmassoption
            disp.solutesoption = self.solutesoption
            # This is in the meantime that ics are read from 
            # what was originally called the dispersion file
            disp.initialconditions = self.initialconditions
            disp.nics = self.nics
            disp.fluxconditions = self.fluxconditions
            disp.npfs = self.npfs

        # Extract model and assign default filenames
        model = args[0]
        if dispersionfilename is None:
            # .rw extension makes more sense 
            # as in fact this file contains most of the 
            # config for rwpt, besides dispersion
            dispersionfilename = f"{model.name}.dispersion"
        self.dispersionfilename = dispersionfilename
        if reconstruction:
            if reconstructionfilename is None:
                reconstructionfilename = f"{model.name}.gpkde"
            self.reconstructionfilename = reconstructionfilename
        self.reconstruction = reconstruction

        # Observations
        if observations is not None:
            if isinstance(observations,list):
                for pobs in observations:
                    if not isinstance(pobs, ModpathRwptObs):
                        raise TypeError('Object in list is type ', type(pobs), '. Expected ModpathRwptObs.')
                # Survived so continue
                self.nobs = len(observations)
                self.observations = observations
            elif isinstance( observations, ModpathRwptObs ):
                self.nobs = 1
                self.observations = [observations]
            else:
                raise TypeError('Observations argument should be of type list or ModpathRwptObs. ', type(observations), ' given.')
        else:
            self.nobs = 0
            self.observations = None


        # done!


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

        Note: mostly the same as classic modpath but with management of rwpt options
        """
    
        f = open(self.fn_path, "w")
        # item 0
        f.write(f"{self.heading}\n")
        # item 1
        f.write(f"{self.mp_name_file}\n")
        # item 2
        f.write(f"{self.listingfilename}\n")
        # item 3
        f.write(
            "{} {} {} {} {} {} {} {} {}\n".format(
                self.simulationtype,
                self.trackingdirection,
                self.weaksinkoption,
                self.weaksourceoption,
                self.budgetoutputoption,
                self.tracemode,
                self.timeseriesoutputoption,
                self.particlesmassoption,
                self.solutesoption,
            )
        )
        # item 4
        f.write(f"{self.endpointfilename}\n")
        # item 5
        if self.simulationtype == 2 or self.simulationtype == 4 or self.simulationtype == 6 :
            f.write(f"{self.pathlinefilename}\n")
        # item 6
        if self.simulationtype == 3 or self.simulationtype == 4 or self.simulationtype == 5 or self.simulationtype == 6:
            f.write(f"{self.timeseriesfilename}\n")
        # item 7 and 8
        if self.tracemode == 1:
            f.write(f"{self.tracefilename}\n")
            f.write(
                f"{self.traceparticlegroup + 1} {self.traceparticleid + 1}\n"
            )
        # item 9
        f.write(f"{self.BudgetCellCount}\n")
        # item 10
        if self.BudgetCellCount > 0:
            v = Util2d(
                self.parent,
                (self.BudgetCellCount,),
                np.int32,
                self.budgetcellnumbers.array + 1,
                name="temp",
                locat=self.unit_number[0],
            )
            f.write(v.string)
    
        # item 11
        f.write(f"{self.referencetimeOption}\n")
        if self.referencetimeOption == 1:
            # item 12
            f.write(f"{self.referencetime[0]:g}\n")
        elif self.referencetimeOption == 2:
            # item 13
            f.write(
                "{:d} {:d} {:g}\n".format(
                    self.referencetime[0] + 1,
                    self.referencetime[1] + 1,
                    self.referencetime[2],
                )
            )
        # item 14
        f.write(f"{self.stoptimeoption}\n")
        if self.stoptimeoption == 3:
            # item 15
            f.write(f"{self.stoptime:g}\n")
    
        # item 16
        if (
                self.simulationtype == 3 or
                self.simulationtype == 4 or
                self.simulationtype == 5 or
                self.simulationtype == 6 ):
            f.write(f"{self.timepointoption}\n")
            if self.timepointoption == 1:
                # item 17
                f.write(
                    f"{self.timepointdata[0]} {self.timepointdata[1][0]}\n"
                )
            elif self.timepointoption == 2:
                # item 18
                f.write(f"{self.timepointdata[0]}\n")
                # item 19
                tp = self.timepointdata[1]
                v = Util2d(
                    self.parent,
                    (tp.shape[0],),
                    np.float32,
                    tp,
                    name="temp",
                    locat=self.unit_number[0],
                )
                f.write(v.string)
    
        # item 20
        f.write(f"{self.zonedataoption}\n")
        if self.zonedataoption == 2:
            # item 21
            f.write(f"{self.stopzone}\n")
            # item 22
            f.write(self.zones.get_file_entry())
    
        # item 23
        f.write(f"{self.retardationfactoroption}\n")
        if self.retardationfactoroption == 2:
            # item 24
            f.write(self.retardation.get_file_entry())
    
        # item 25
        if self.particlegroups is None:
            f.write(f"0\n")
        else:
            f.write(f"{len(self.particlegroups)}\n")
            if ( self.particlesmassoption == 0):
                for pg in self.particlegroups:
                    pg.write(f, ws=self.parent.model_ws)
            elif ( self.particlesmassoption == 1):
                for pg in self.particlegroups:
                    pg.write(f, ws=self.parent.model_ws, mass=True)
            elif ( self.particlesmassoption == 2):
                for pg in self.particlegroups:
                    pg.write(f, ws=self.parent.model_ws, mass=True, solute=True)
   

        # MODPATH-RW
        if self.simulationtype > 4:

            # RWPT config filename
            f.write(f"{self.dispersionfilename}\n")

            if self.reconstruction: 
                # Possibly requires timeseries 
                #if self.simulationtype == 7:
                #    raise Exception('')
                f.write(f"1\n")
                f.write(f"{self.reconstructionfilename}\n")
            else:
                f.write(f"0\n")
      

            if self.nobs > 0:
                f.write(f"{len(self.observations)}\n")
                # Write each
                for obs in self.observations:
                    obs.write(f=f)
            else:
                pass


        # And close
        f.close()




# Needs transformation of solute concentrations into particles

