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
        longitudinal = 1,
        transverse   = 0.1,
        molecular    = 0.0,
        advection    = 'eulerian',
        dimensions   = 3,
        timestep     = 'courant',
        courant      = 0.1,
        peclet       = 100.0,
        ctdisp       = 0.1,
        extension    = 'dispersion',
        icbound      = 0,
        initialconditions = None, 
        mass              = 1.0,
        icformat          = 2, # 1: density, 2: classic particles
        massformat        = 1, # 1:
        injectioncells    = None, 
        injectionseries   = None, 
        injectionmass     = 1.0,
        modelkind         = None, # linear, nonlinear
        betatransverse    = 0.5, 
        betalongitudinal  = 1,
        mediumdelta       = 5, 
        mediumdsize       = 1,
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
        #self.parent._generate_heading()


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
        self.icbound= Util3d(
            model,
            shape3d,
            np.int32,
            icbound,
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
        self.ic = initialconditions


        # EXPECTED TO HAVE THE SAME LENGTH
        if self.ic is not None:
            nics = len(self.ic)
            self.mass     = mass 
            self.icformat = icformat
            self.massformat = massformat
            if isinstance(self.mass, (float,int)):
                self.mass = np.repeat( self.mass, nics )
            if isinstance(self.icformat, (float,int)):
                self.icformat = np.repeat( self.icformat, nics ) 
            if isinstance(self.massformat, (float,int)):
                self.massformat = np.repeat( self.massformat, nics ) 


        # INJECTION BOUNDARY CONDITIONS
        self.injectioncells = injectioncells
        self.injectionseries = injectionseries
        
        if injectioncells is not None:
            self.injectionformat = []
            for icell in injectioncells:
                if (
                        (injectionseries is None) or 
                        ( len(injectionseries) != len(injectioncells) )
                    ):
                    self.injectionformat.append(1) 
                else:
                    # TIMESERIES
                    self.injectionformat.append(2) 

            nijs = len( injectioncells ) 
            self.injectionmass = injectionmass 
            if isinstance(self.injectionmass, (float,int)):
                self.injectionmass = np.repeat( self.injectionmass, nijs )




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
        #f.write(f"# {self.heading}\n")

        
        # Select model dispersivity
        if self.modelkind == 1:

            # Write modelkind
            f.write(f"{self.modelkind}\n")

            # Write dispersivities
            f.write(self.longitudinal.get_file_entry())
            f.write(self.transverse.get_file_entry())
            # Write molecular diffusion
            f.write(f"{self.molecular:.16f}\n")

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

        # Write icbound
        f.write(self.icbound.get_file_entry())

        # Write initial conditions
        if self.ic is not None:
            f.write(f"{len(self.ic)}\n")
            for idpg, pg in enumerate(self.ic):
                f.write(f"{self.icformat[idpg]}\n")

                if self.icformat[idpg] == 1:
                    f.write(f"MPG{idpg}\n")
                    f.write(f"{self.massformat[idpg]}\n")

                    if self.massformat[idpg] == 1:
                        f.write(f"{self.mass[idpg]:.10f}\n")      # pggroup should carry this
                        icarray = Util3d(
                            self.model,
                            self.shape3d,
                            np.float32,
                            pg,
                            name="IC"+str(idpg),
                            locat=self.unit_number[0],
                        )

                        # Write icarray
                        f.write(icarray.get_file_entry())
                    else:
                        pass

                elif self.icformat[idpg] == 2:
                    pg.write(f, ws=self.parent.model_ws)
                    f.write(f"{self.mass[idpg]:.10f}\n")      # pggroup should carry this
        else:
            f.write(f"0\n")


        if self.injectioncells is not None:
            # Write injection cells information
            f.write(f"{len(self.injectioncells)}\n")
            for idpg, pg in enumerate(self.injectioncells):
                f.write(f"{self.injectionformat[idpg]}\n")    # FORMAT
                f.write(f"IPG{idpg}\n")                       # NAME
                f.write(f"{self.injectioncells[idpg][0]}\n")  # CELLNUMBER
                f.write(f"{self.injectionmass[idpg]:.10f}\n") # INJECTIONMASS
                if self.injectionformat[idpg] == 1:
                    f.write(f"{self.injectioncells[idpg][1]}\n")  # CONCENTRATION
                    # SOMETHING LIKE A FINAL TIME OR SOMETHING
                    # INITiALTIME AND FINALTIME
                elif self.injectionformat[idpg] == 2: # AS TIMESERIES
                    f.write(f"{len(self.injectionseries[idpg][0])}\n")  # CONCENTRATION
                    df = pd.DataFrame()
                    df['time'] = self.injectionseries[idpg][0] 
                    df['conc'] = self.injectionseries[idpg][1]
                    filename   = 'injcell'+str(idpg)+'_timeseries.csv'
                    df.to_csv(os.path.join(
                        self.model.model_ws, filename),
                        header=None,
                        index=False,
                        sep=' ', 
                        float_format="%.10f"
                    )
                    f.write(f"{filename}\n")
        else:
            f.write(f"0\n")

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
        skiptimeserieswriter = True,
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
                f.write(f"{b:10f}")

        # Domain sizes
        for idb, b in enumerate(self.domainsize):
            if idb == len(self.domainsize)-1: 
                f.write(f"{b:10f}\n")
            else:
                f.write(f"{b:10f}")

        # Bin sizes
        for idb, b in enumerate(self.binsize):
            if idb == len(self.binsize)-1: 
                f.write(f"{b:10f}\n")
            else:
                f.write(f"{b:10f}")

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
                    f.write(f"{b:10f}")
        else:
            # Brute force reconstruction
            f.write(f"0\n") # 0 for id into fortran
            
            # kernel params: minh/lambda, maxh/lambda
            for idb, b in enumerate(self.kerneldatabaseparams):
                if idb == 0:
                    f.write(f"{b:16f}")
                if idb == 2: 
                    f.write(f"{b:16f}\n")


        # And close
        f.close()



#class ModpathRwptInitialCondition( Package ):
#    """
#    MODPATH RWPT InitialCondition Package Class.
#    Parameters
#    ----------
#    model : model object
#        The model object (of type :class:`flopy.modpath.Modpath7`) to which
#        this package will be added.
#    binsize : list [sx,sy,sz]
#        List with cell size employed for reconstruction process 
#    domainsize : list [sx,sy,sz]
#        List with domain sizes. These are used for computing reconstruction 
#        grid dimensions using given binsize.
#    noptloops: int
#        The number of optimization loops for determinin smoothing conditions
#    outputfilename: str
#        Filename that will be used for writing reconstruction output. Contains 
#        cell id, density and both time step and particle group identification.
#    extension : str, optional
#        File extension (default is 'dispersion').
#    """
#
#    def __init__(
#        self,
#        model,
#        binsize        = [1,1,1],
#        domainsize     = [1,1,1],
#        noptloops      = 10, 
#        outputfilename = 'gpkde.output',
#        extension      = 'gpkde',
#    ):
#
#        unitnumber = model.next_unit()
#
#        super().__init__(model, extension, "IC", unitnumber)
#
#        shape = model.shape
#        if len(shape) == 3:
#            shape3d = shape
#        elif len(shape) == 2:
#            shape3d = (shape[0], 1, shape[1])
#        else:
#            shape3d = (1, 1, shape[0])
#
#        # Parameters assignment
#        self.binsize    = binsize
#        self.domainsize = domainsize
#        self.noptloops  = noptloops
#
#        if outputfilename is not None:
#            self.outputfilename = outputfilename
#        else:
#            raise Exception('requires outputfilename for reconstruction output')
#
#        self.parent.add_package(self)
#
#        
#
#    def write_file(self, check=False):
#        """
#        Write the package file
#        Parameters
#        ----------
#        check : boolean
#            Check package data for common errors. (default False)
#        Returns
#        -------
#        None
#        """
#        # Open file for writing
#        f = open(self.fn_path, "w")
#
#        # Output filename
#        f.write(f"{self.outputfilename}\n")
#
#        # Domain sizes
#        for idb, b in enumerate(self.domainsize):
#            if idb == len(self.domainsize)-1: 
#                f.write(f"{b:16f}\n")
#            else:
#                f.write(f"{b:16f}")
#
#        # Bin sizes
#        for idb, b in enumerate(self.binsize):
#            if idb == len(self.binsize)-1: 
#                f.write(f"{b:16f}\n")
#            else:
#                f.write(f"{b:16f}")
#
#
#        # Number of optimization loops 
#        f.write(f"{self.noptloops:10d}\n")
#
#        # more to come...
#
#        # And close
#        f.close()








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
            reconstruction=False,
            dispersionfilename=None,
            reconstructionfilename=None,
            **kwargs
        ):

        # Call parent constructor
        super().__init__(*args,**kwargs, extension='mprwpt') 

        # Extract model and assign default filenames
        model = args[0]
        if dispersionfilename is None:
            dispersionfilename = f"{model.name}.dispersion"
        self.dispersionfilename = dispersionfilename
        if reconstruction:
            if reconstructionfilename is None:
                reconstructionfilename = f"{model.name}.gpkde"
            self.reconstructionfilename = reconstructionfilename
        self.reconstruction = reconstruction

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

        DEV: same as upstream file but with specific management of rwpt options
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
            "{} {} {} {} {} {}\n".format(
                self.simulationtype,
                self.trackingdirection,
                self.weaksinkoption,
                self.weaksourceoption,
                self.budgetoutputoption,
                self.tracemode,
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
            # ORIGINAL
            #f.write(f"{self.stoptime + 1:g}\n")
    
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
        f.write(f"{len(self.particlegroups)}\n")
        for pg in self.particlegroups:
            pg.write(f, ws=self.parent.model_ws)
   

        # RWPT
        if self.simulationtype > 4:
            f.write(f"{self.dispersionfilename}\n")

            if self.reconstruction: 
                # Possibly requires timeseries 
                #if self.simulationtype == 7:
                #    raise Exception('')
                f.write(f"1\n")
                f.write(f"{self.reconstructionfilename}\n")
            else:
                f.write(f"0\n")
       

        # And close
        f.close()




# Needs transformation of solute concentrations into particles

