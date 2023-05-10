'''
Implement classes and methods for configuring MODPATH-RW sims with flopy
'''

# Python
import os
import numpy as np
import pandas as pd
from enum import Enum
import warnings

# Flopy
import flopy
import flopy.modpath.mp7sim as mp7
from flopy.pakbase import Package
from flopy.utils import Util2d, Util3d
from flopy.modpath.mp7particlegroup import (
    ParticleGroup,
    ParticleGroupLRCTemplate,
    ParticleGroupNodeTemplate,
)
from flopy.modpath.mp7particledata import (
    ParticleData,
)

from .utils import count_instances


# Updates enumeration of available simulations
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


# Update particle groups classes
# Extend interfaces to understand mass and solute
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


class ModpathRW( flopy.modpath.Modpath7 ):
    '''
    MODPATH-RW class 

    Extends flopy.modpath.Modpath7
    Overloads write_name_file
    '''
    def __init__(
            self, *args,
            version='modpathrw',
            simfile_ext="mprw" , # Not being recognized
            **kwargs
        ):

        # Call parent constructor
        super().__init__(*args,**kwargs)
        self.version_types = {"modpathrw": "MODPATH-RW"}
        self.set_version(version)

        # Following filenames are generated
        # after parent constructor, appended to the class
        #   - mpbas_file
        #   - dis_file
        #   - grbdis_file
        #   - tdis_file
        #   - headfilename
        #   - budgetfilename

        # Specific MODPATH-RW filenames are 
        # extracted from their respective package ( if given )


    # Overload repr
    def __repr__(self):
        return "MODPATH-RW model"


    # Overload the write_name_file method
    def write_name_file(self):
        """
        Write the name file

        Returns
        -------
        None

        """
        fpth = os.path.join(self.model_ws, self.mpnamefile)
        f = open(fpth, "w")
        f.write(f"{self.heading}\n")
        if self.mpbas_file is not None:
            f.write(f"MPBAS      {self.mpbas_file}\n")
        if self.dis_file is not None:
            f.write(f"DIS        {self.dis_file}\n")
        if self.grbdis_file is not None:
            f.write(f"{self.grbtag:10s} {self.grbdis_file}\n")
        if self.tdis_file is not None:
            f.write(f"TDIS       {self.tdis_file}\n")
        if self.headfilename is not None:
            f.write(f"HEAD       {self.headfilename}\n")
        if self.budgetfilename is not None:
            f.write(f"BUDGET     {self.budgetfilename}\n")

        # MODPATH-RW specifc files
        dsppkg = self.get_package('DSP')
        if dsppkg is not None:
            f.write(f"DSP        {dsppkg.file_name[0]}\n")
        spcpkg = self.get_package('SPC')
        if spcpkg is not None:
            f.write(f"SPC        {spcpkg.file_name[0]}\n")
        rwoptspkg = self.get_package('RWOPTS')
        if rwoptspkg is not None:
            f.write(f"RWOPTS     {rwoptspkg.file_name[0]}\n")
        gpkdepkg = self.get_package('GPKDE')
        if gpkdepkg is not None:
            f.write(f"GPKDE      {gpkdepkg.file_name[0]}\n")
        icpkg = self.get_package('IC')
        if icpkg is not None:
            f.write(f"IC         {icpkg.file_name[0]}\n")
        bcpkg = self.get_package('BC')
        if bcpkg is not None:
            f.write(f"BC         {bcpkg.file_name[0]}\n")
        srcpkg = self.get_package('SRC')
        if srcpkg is not None:
            f.write(f"SRC        {srcpkg.file_name[0]}\n")
        obspkg = self.get_package('OBS')
        if obspkg is not None:
            f.write(f"OBS        {obspkg.file_name[0]}\n")

        # Done
        f.close()

        return


class ModpathRWSim( flopy.modpath.Modpath7Sim ):
    '''
    MODPATH-RW Simulation File Package Class. 

    Extends from Modpath7Sim 
    '''

    def __init__( self, *args,
            timeseriesoutputoption=0,
            particlesmassoption=0,
            solutesoption=0,
            fluxconditions=None, # DEPRECATE
            **kwargs
        ):

        # Call parent constructor
        super().__init__(*args,**kwargs, extension='mprw' )


        # Extract model
        model = args[0]


        # Override interpretation of particlegroups, allows None
        if 'particlegroups' not in kwargs.keys():
            particlegroups = None
            self.particlegroups = particlegroups


        # New options
        if (timeseriesoutputoption not in [0,1,2]):
            raise ValueError('Timeseries output option should be 0, 1 or 2. Given :', str(timeseriesoutputoption) )
        self.timeseriesoutputoption = timeseriesoutputoption
        if (particlesmassoption not in [0,1,2]):
            raise ValueError('Particles mass option should be 0, 1 or 2. Given :', str(particlesmassoption) )
        self.particlesmassoption = particlesmassoption
        if (solutesoption not in [0,1]):
            raise ValueError('Solutes option should be 0 or 1. Given :', str(solutesoption) )
        self.solutesoption = solutesoption


        # NOT IMPLEMENTED
        # Prescribed flux conditions
        if fluxconditions is not None:
            if isinstance(fluxconditions,list):
                for pic in fluxconditions:
                    if not isinstance(pic, ModpathRWFlux):
                        raise TypeError('Object in list is type ', type(pic), '. Expected ModpathRWFlux.')
                # Survived so continue
                self.npfs = len(fluxconditions)
                self.fluxconditions = fluxconditions
            elif isinstance( fluxconditions, ModpathRWFlux ):
                self.npfs = 1
                self.fluxconditions = [fluxconditions]
            else:
                raise TypeError('Initial conditions argument should be of type list or ModpathRWFlux. ', type(fluxconditions), ' given.')
        else:
            self.npfs = 0
            self.fluxconditions = None


        # If simulation is RW
        if self.simulationtype > 4: 

            # Assign some properties to parent obj
            # Needed by: SPC, IC
            self._parent.particlesmassoption = particlesmassoption
            self._parent.solutesoption = solutesoption


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
   
        # And close
        f.close()

        return



class ModpathRWBas( flopy.modpath.Modpath7Bas ):
    '''
    MODPATH-RW Basic class 
    Extends flopy.modpath.Modpath7Bas
    '''
    def __init__(self, *args, **kwargs):
        # Call parent constructor
        super().__init__(*args,**kwargs)



class ModpathRWOptions( Package ):
    """
    MODPATH-RW Options Package Class.
    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    timestep : str
        Define method for computing particles time step. Can be adv, disp, min or fixed. 
    courant  : float
        Courant number that will be used at each cell for computing time step with adv criteria 
    ctdisp   :  float 
        CT constant for computing time step with disp criteria
    deltat   :  float 
        Value employed when timestep selection is fixed
    advection : str
        Advection model ot be used for advective component in RWPT integration.
        Could be set to 'eulerian' or 'exponential' (default is 'eulerian').
    dimensionsmask : [int,int,int]
        Determine on which dimensions RW displacements are calculated. 
        A value of 1 indicates that dimension is active and 0 inactive.
        For example, for a 2D RW model (x,y model and single layer), dimensionsmask should be [1,1,0].
        At least one dimension is required. Defaults to 3D ( [1,1,1] )
    extension : str, optional
        File extension (default is 'dispersion').
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

        # model shape ( needed ? )
        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])
        self.model = model
        self.shape3d = shape3d

        # timestep
        if ( timestep not in ['adv', 'disp', 'min', 'fixed'] ):
            raise Exception('flopyrw:ModpathRWOptions: time step selection model ' + timestep + ' is not valid. Should be: adv, disp, min, fixed.')
        self.timestep  = timestep
        # Trust
        self.courant   = courant
        self.ctdisp    = ctdisp 
        self.deltat    = deltat

        # advection 
        if ( advection not in ['eulerian','exponential'] ):
            raise Exception('flopyrw:ModpathRWOptions: advection model ' + advection + ' is not valid. Use eulerian or exponential')
        self.advection  = advection

        # dimensionsmask
        if not isinstance( dimensionsmask, (list, np.array) ):
            raise Exception('flopyrw:ModpathRWOptions: dimensionsmask should be list or np.array. Given ', type(dimensionsmask) )
        if isinstance( dimensionsmask, list ):
            if len(dimensionsmask) != 3:
                raise Exception('flopyrw:ModpathRWOptions: dimensionsmask list should be of len(dimensionsmask) = 3. Given ', len(dimensionsmask) )
            dimensionsmask = np.array(dimensionsmask).astype(np.int32)
        if len(dimensionsmask) != 3:
            raise Exception('flopyrw:ModpathRWOptions: dimensionsmask should be of len(dimensionsmask) = 3. Given ', len(dimensionsmask) )
        for nd in range(3):
            if (dimensionsmask[nd] not in [0,1] ):
                raise Exception('flopyrw:ModpathRWOptions: values in dimensionsmask should be 0 or 1. At position ', str(nd), ' given ', str(dimensionsmask[nd]) )
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
        else:
            raise Exception('flopyrw:ModpathRWOptions: time step selection model ' + timestep + ' is not valid. Should be: adv, disp, min, fixed.')


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
        


class ModpathRWGpkde( Package ):
    """
    MODPATH-RW GPKDE Package Class.
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
        The number of optimization loops smoothed reconstruction
    asconcentration: bool
        Flag to indicate whether reconstruction should be returned as resident concentration. 
        Because the reconstruction grid is different than flowmodel grid, said transformation 
        can be easily achieved in case both porosities and retardation are spatially uniform.
        If not, then some kind of intersection or relation should be established between the 
        base flowmodel and the reconstruction grid.
    outputfilename: str
        Filename that will be used for writing reconstruction output. Contains 
        cell id, density and both time step and particle group identification.
    extension : str, optional
        File extension (default is 'dispersion').
    """

    def __init__(
        self,
        model,
        binsize         = [1,1,1],
        domainsize      = [1,1,1],
        domainorigin    = None, 
        minhlambda      = 1.0 ,
        maxhlambda      = 0.1 ,
        deltahlambda    = 10.0,
        kerneldatabase  = False, 
        noptloops       = 2,
        asconcentration = False,
        outputfilename  = None,
        extension       = 'gpkde',
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

        # Some health checks
        self.binsize              = np.array(binsize).astype(np.float32)
        self.domainsize           = np.array(domainsize).astype(np.float32)

        # Fills domain origin with some values 
        # detemined from the modelgrid definition
        if domainorigin is None :
            domainorigin = [
                    self._parent.flowmodel.modelgrid.xoffset,
                    self._parent.flowmodel.modelgrid.yoffset,
                    np.min(self._parent.flowmodel.modelgrid.botm)
                ]
        self.domainorigin         = np.array(domainorigin).astype(np.float32)


        # Depending on the kind of modelgrid maybe some
        # decisions can be made for binsizes and/or domain sizes
        # StructuredGrid has the is_regular function

        self.noptloops            = noptloops
        self.kerneldatabase       = kerneldatabase       
        self.minhlambda           = minhlambda
        self.maxhlambda           = maxhlambda
        self.deltahlambda         = deltahlambda

        self.asconcentration      = asconcentration

        if outputfilename is not None:
            self.outputfilename = outputfilename
        else:
            self.outputfilename = f"{model.name}.gpkde.out"

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

        # Open file for writing
        f = open(self.fn_path, "w")

        # Output filename
        f.write(f"{self.outputfilename}\n")

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

        # Reconstruction with kernel database or brute force
        # 1: kerneldatabase
        # 0: bruteforcedensity
        if self.kerneldatabase:
            # Kernel database reconstruction
            f.write(f"1\n") # 1 for id into fortran

            # Database params: minh/lambda, daltah/lambda, maxh/lambda
            f.write(f"{self.minhlambda:10f}\t")
            f.write(f"{self.deltahlambda:10f}\t")
            f.write(f"{self.maxhlambda:10f}\n")

        else:
            # Brute force reconstruction
            f.write(f"0\n") # 0 for id into fortran
            
            # kernel params: minh/lambda, maxh/lambda
            f.write(f"{self.minhlambda:10f}\t")
            f.write(f"{self.maxhlambda:10f}\n")


        if self.asconcentration:
            f.write(f"1\n") # 1 for id into fortran
        else:
            f.write(f"0\n") # 0 for id into fortran


        # And close
        f.close()

        return





class ModpathRWSpc( Package ):
    """
    MODPATH-RW Species Package Class.
    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    dispersion: ModpathRWDsp object
        Defines dispersion model and properties to this specie.
    pgroups: list, np.array
        Particle groups related to this specie.
    """


    # Class properties
    COUNTER    = 0
    UNITNUMBER = 0
    INSTANCES  = []


    @count_instances
    def __init__(
        self,
        model,
        dispersion, # it's like a foreign key
        pgroups    = None,
        id         = None,
        stringid   = None,
        extension  = 'spc',
    ):
        
        # Define UNITNUMBER if the first instance is created
        if self.__class__.COUNTER == 1:
            unitnumber = model.next_unit()
            super().__init__(model, extension, "SPC", unitnumber)
            self.__class__.UNITNUMBER = self.unit_number[0]
        else:
            # Needed for Util3d
            self._parent = self.INSTANCES[0]._parent


        ## Necessary ?
        #shape = model.shape
        #if len(shape) == 3:
        #    shape3d = shape
        #elif len(shape) == 2:
        #    shape3d = (shape[0], 1, shape[1])
        #else:
        #    shape3d = (1, 1, shape[0])
        #self.model = model
        #self.shape3d = shape3d


        # Assign dispersion
        if not isinstance( dispersion, ModpathRWDsp ):
            raise Exception('flopyrw:ModpathRWSpc: dispersion parameter should be of ModpathRWDsp type. Given ', type( dispersion ) )
        self.dispersion = dispersion

        # Assign pgroups id's
        if pgroups is None: 
            warnings.warn("flopyrw:ModpathRWSpc: no pgroups were specified, defaults to pgroups = [0]")
            pgroups = [0]
            # Or leave as None
        if not isinstance(pgroups, (list,np.array)):
            raise Exception('flopyrw:ModpathRWSpc: pgroups should be of type list or np.array. Given ', type( pgroups ) )
        if isinstance( pgroups, list ):
            pgroups = np.array( pgroups )
        # Requires some more sanity checks to verify integer id's and so on
        # Might be worth considering removing particle groups altogether and extend 
        # species from those classes.
        self.pgroups = pgroups
        

        # Define id
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER

        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = 'SPC'+str(self.__class__.COUNTER)


        # Add package and save instance        
        if self.__class__.COUNTER == 1:
            self.parent.add_package(self)
        self.__class__.INSTANCES.append( self )


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

        # Open file for writing
        f = open(self.INSTANCES[0].fn_path, "w")

        # Write how many INSTANCES 
        # and loop over them
        f.write(f"{self.COUNTER}\n")

        for ins in self.INSTANCES:

            # Write id
            f.write(f"{ins.id}\n")

            # Write stringid
            f.write(f"{ins.stringid}\n")

            # Need to write the related groups 
            # only if particlesmassoption not equal 2
            # If particlesmassoption equal 2, then 
            # soluteId is infered from given particlegroups
            if ( self.INSTANCES[0]._parent.particlesmassoption != 2 ):
                # Write the number of groups
                f.write(f"{len(ins.pgroups)}\n")

                for idpg, pg in enumerate(ins.pgroups):
                    # Write the pgroup id in the list of pgroups
                    f.write(f"{pg+1}\n")

            # Write dispersion "foreign key"        
            f.write(f"{ins.dispersion.id}\n")
            f.write(f"{ins.dispersion.stringid}\n")


        # Done
        return


class ModpathRWIc( Package ):
    """
    MODPATH-RW Initial Condition Package Class.
    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    concentration : ndarray
        Interpreted as resident concentration, that is, the product beween porosity and
        dissolved concentration. 
    """


    # Class properties
    COUNTER    = 0
    UNITNUMBER = 0
    INSTANCES  = []


    @count_instances
    def __init__(
        self,
        model,
        kind          = 1,
        mass          = .1,
        concentration = None,
        soluteid      = 1, # replace by specie ?
        id            = None,
        stringid      = None,
        extension     = 'ic',
    ):
       

        # Define UNITNUMBER if the first instance is created
        if self.__class__.COUNTER == 1:
            unitnumber = model.next_unit()
            super().__init__(model, extension, "IC", unitnumber)
            self.__class__.UNITNUMBER = self.unit_number[0]
        else:
            # Needed for Util3d
            self._parent = self.INSTANCES[0]._parent


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
            raise ValueError('flopyrw:ModpathRWIc: Invalid value for kind, should be 1. Given ', str(kind) )
        self.kind = kind

        # Need health checks
        self.mass       = mass
        self.soluteid   = soluteid

        # define id
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER

        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = 'IC'+str(self.__class__.COUNTER)

        # and concentration: needs id
        if (
            ( concentration is not None ) and 
            ( self.kind == 1 )
        ):
            if not isinstance( concentration, np.ndarray ):
                raise TypeError('flopyrw:ModpathRWIc: concentration should be an np.ndarray')
            self.concentration = Util3d(
                model,
                shape3d,
                np.float32,
                concentration,
                name=self.stringid,
                locat=self.__class__.UNITNUMBER,
            )
        elif ( (self.kind ==1) and (concentration is None) ):
            raise ValueError('flopyrw:ModpathRWIc: kind 1 requires concentration array, None given.')
     

        # Add package and save instance        
        if self.__class__.COUNTER == 1:
            self.parent.add_package(self)
        self.__class__.INSTANCES.append( self )


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


        # Open file for writing
        f = open(self.INSTANCES[0].fn_path, "w")

        # Write how many INSTANCES 
        # and loop over them
        f.write(f"{self.COUNTER}\n")

        for ins in self.INSTANCES:

            # Write string id 
            f.write(f"{ins.stringid}\n")

            # Kind/format of initial condition 
            f.write(f"{ins.kind}\n")

            # 1: resident concentration array  
            if ins.kind == 1:

                # Give particles mass 
                f.write(f"{ins.mass:.10f}\n")

                # If solutes are being specified, do it
                if( 
                    ( self.INSTANCES[0]._parent.particlesmassoption == 2 ) or
                    ( self.INSTANCES[0]._parent.solutesoption == 1 ) 
                ):
                    f.write(f"{ins.soluteid}\n")

                # And concentration distribution
                f.write(ins.concentration.get_file_entry())


        # Done
        return






from collections import Counter
class ModpathRWSrc( Package ):
    """
    MODPATH-RW Source Package Class.
    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    """


    # Class properties
    COUNTER = 0
    UNITNUMBER = 0
    INSTANCES = []


    @count_instances
    def __init__(
        self,
        model,
        format        = 'aux',
        sources       = None,
        cells         = None, 
        mass          = 1.0,
        nparticles    = 2,
        concentration = None, 
        soluteid      = None,
        sourcename    = None, 
        auxiliary     = None,
        id            = None, 
        stringid      = None, 
        extension     = 'src',
    ):
       
        # Define UNITNUMBER if the first instance is created
        if self.__class__.COUNTER == 1:
            unitnumber = model.next_unit()
            super().__init__(model, extension, "SRC", unitnumber)
            self.__class__.UNITNUMBER = self.unit_number[0]
        else:
            # Needed for Util3d
            self._parent = self.INSTANCES[0]._parent

        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])
        self.model = model
        self.shape3d = shape3d

     
        self.format = format
        self.modelversion = self._parent.flowmodel.version 
       

        if mass is not None:
            DEFAULTMASS = mass
        else: 
            DEFAULTMASS = 1.0

        if nparticles is not None:
            DEFAULTNP = nparticles
        else:
            DEFAULTNP = 2

        # Clarify whether this is python based or fortran based
        if soluteid is not None:
            DEFAULTSOLID = int(soluteid)
        else:
            DEFAULTSOLID = int(1)


        # If format is aux, validate sources, verify auxiliary variables 
        sourceslist = [] # STORE IN CLASS ?
        alldatalist = []

        # Validate format kind
        if ( ( format.upper() not in ['AUX','AUXILIARY'] ) ):
            raise Exception('flopyrw:ModpathRWSrc: source format ' + format.upper() + ' not implemented. ')

        # Validate that sources was given
        if ( ( sources is None ) and ( format.upper() in ['AUX','AUXILIARY'] ) ):
            raise Exception('flopyrw:ModpathRWSrc: package requires sources data. None was given.')

        # Read sources
        for sd in sources:

            pname = sd[0]
            if not isinstance( pname, str ): 
                raise TypeError('flopyrw:ModpathRWSrc: package name should be str. ',type(pname),' was given.')
            pkg = self._parent.flowmodel.get_package(pname.lower())
            if pkg is None:
                raise Exception(\
                    'flopyrw:ModpathRWSrc: package name '+ pname.lower() \
                    + ' was not found in flowmodel.'+' Available packages are: ' +",".join(self._parent.flowmodel.get_package_list()) )

            # Get the auxlist from package
            # Note: auxlist filled only with str.upper()

            if self.modelversion == 'mf6':
                # Get the aux data and create auxlist
                if not pkg.auxiliary.has_data():
                    raise Exception('flopyrw:ModpathRWSrc: auxiliary property does not has data for package ' + pkg.name[0])
                aux     = pkg.auxiliary.array[0]
                auxlist = [a.upper() for a in aux ] 
            else:
                # Extract possible aux variable names and perform some preliminary checks for != MF6
                opts = pkg.options
                if isinstance(opts, flopy.utils.OptionBlock):
                    auxlistcand = opts.auxillary
                elif isinstance(opts,(list,tuple)):
                    auxlistcand = pkg.options
                else:
                    raise TypeError('flopyrw:ModpathRWSrc: unexpected options type in package ' + pkg.name[0] )

                if ( len(auxlistcand) == 0):
                    raise Exception('flopyrw:ModpathRWSrc: no auxiliary variables were found in package ' + pkg.name[0] )

                # Loop over possible aux vars and create auxlist
                auxlist = []
                for ia, auxc in enumerate(auxlistcand):
                    # Clean the name ( str is assumed )

                    # Clean borders
                    cand = auxc.strip().upper()
                    
                    # If is the identifier str, go the next
                    if ( cand in ['AUX','AUXILIARY'] ):
                        # To the next
                        continue

                    # Break it
                    cand = cand.split()

                    # It it had a whitespace
                    if ( len(cand) == 2 ):
                        # Need to verify that the first word 
                        # is the aux identifier and the second
                        # is the auxvar name
                        if ( cand[0].upper() in ['AUX','AUXILIARY'] ):
                            # Alles gud
                            auxlist.append( cand[1].upper() )
                    elif ( len(cand) == 1):
                        auxlist.append( cand[0].upper() )
                    else:
                        # Probably the case in which auxvars have whitespace in between. omg
                        print( cand )
                        raise Exception('flopyrw:ModpathRWSrc: something wrong with definition of aux in package ' + pkg.name[0])

            # Determine the format in which the aux variables are specified
            auxspec = sd[1]
            if isinstance( auxspec, str ):
                # If a str was given, transform to list
                auxnames     = [auxspec.upper()]
                auxmasses    = [DEFAULTMASS]
                auxtemplates = DEFAULTNP*np.ones(shape=(1,3),dtype=np.int32)
                auxsolutes   = [DEFAULTSOLID] # Clarify whether these are python based or fortran based
            elif isinstance( auxspec, (list,tuple)):
                # List to store aux specs
                auxnames     = []
                auxmasses    = DEFAULTMASS*np.ones( shape=(len(auxspec),) )
                auxtemplates = DEFAULTNP*np.ones(shape=(len(auxspec),3),dtype=np.int32)
                auxsolutes   = DEFAULTSOLID*np.ones( shape=(len(auxspec),),dtype=np.int32 ) # Clarify whether these are python based or fortran based
                for iaspc, aspc in enumerate(auxspec):
                    if isinstance( aspc, (list,tuple) ):
                        auxspeclen = len(aspc)
                        if auxspeclen == 1:
                            # Consider as auxvar
                            if not isinstance( aspc, str ): 
                                raise TypeError('flopyrw:ModpathRWSrc: auxiliary name should be str. ',type(aspc),' was given.')
                            # Good
                            auxnames.append( aspc.upper() )
                        elif ( (auxspeclen >= 2 ) and (auxspeclen <= 4 ) ):
                            # First is auxvarname
                            if not isinstance( aspc[0], str ): 
                                raise TypeError('flopyrw:ModpathRWSrc: auxiliary name should be str. ',type(aspc[0]),' was given.')
                            # Good
                            auxnames.append( aspc[0].upper() )
                            
                            # Second is mass 
                            if not isinstance( aspc[1], (int,float) ): 
                                raise TypeError('flopyrw:ModpathRWSrc: particles mass should be int/float. ',type(aspc[1]),' was given.')
                            # Good
                            auxmasses[iaspc] = aspc[1]
                           
                            # Continue reading 
                            if auxspeclen >= 3:
                                # Consider as [auxvar,mass,template]
                                if isinstance(aspc[2],(list,tuple)):
                                    for inp, npa in enumerate(aspc[2]):
                                        auxtemplates[iaspc,inp] = npa
                                elif isinstance( aspc[2], int ):
                                    # Assume uniform template 
                                    auxtemplates[iaspc,:] = aspc[2]
                                else:
                                    raise TypeError('flopyrw:ModpathRWSrc: particles template should be int/list/tuple. ',type(aspc[2]),' was given.')

                                if auxspeclen == 4:
                                    # Consider as [auxvar,mass,template,soluteid]
                                    if not isinstance(aspc[3],int):
                                        raise TypeError('flopyrw:ModpathRWSrc: solute id should be int. ',type(aspc[3]),' was given.')
                                    auxsolutes[iaspc] = aspc[3]
                        else:
                            raise ValueError('flopyrw:ModpathRWSrc: max length of aux specification is 4. ' + str(len(aspc)) + ' was given.')
                    elif isinstance( aspc, str ):
                        # Is a str, give to the list
                        auxnames.append( aspc.upper() )
                    else: 
                        raise TypeError('flopyrw:ModpathRWSrc: auxiliary name should be str. ',type(aspc),' was given.')

            else:
                raise TypeError('flopyrw:ModpathRWSrc: auxiliary name should be str or list/tuple. ',type(auxspec),' was given.')


            # Validate given auxnames
            for ian, an in enumerate(auxnames):
                acount = auxlist.count(an.upper())
                if acount > 1:
                    raise Exception('flopyrw:ModpathRWSrc: auxiliary name '+an.upper()+' was specified more than once in package ' + pname.upper() +'. Check definition.')
                if an.upper() not in auxlist:
                    raise ValueError('flopyrw:ModpathRWSrc: auxiliary name ' +an.upper()+' was not specified in package '+ pname.upper() +'. Check package definition.') 

                # If it got until here, definition is consistent,
                # save the pair [sourcename,auxname]
                tempair  = [pname.upper(),an.upper()]
                sourceslist.append(tempair)
                tempdata = [an.upper(), auxmasses[ian], auxtemplates[ian,0], auxtemplates[ian,1], auxtemplates[ian,2], auxsolutes[ian] ]
                alldatalist.append(tempdata)
        ## end read sources ## 

        # Validate that combinations of [sourcename, auxname] are unique 
        counter       =  Counter(frozenset(sl) for sl in sourceslist )
        counterstatus = [counter[frozenset(sl)] == 1 for sl in sourceslist]

        # Break if some inconsistencies
        if not np.all( counterstatus ):
            raise Exception('flopyrw:ModpathRWSrc: some pairs of [sourcename,auxname] are repeated. Combination should be unique for all source packages.')

        # Combination of sourcenames,auxnames is valid
        # Extract unique source names
        self.uniquesources = np.unique( np.array(sourceslist)[:,0] )
        self.allauxnames   = []
        self.alldataperaux = []
        for ius, us in enumerate( self.uniquesources ):
            auxnamespersource = []
            dataperaux        = []
            for isl, sl in enumerate(sourceslist):
                if sl[0] == us:
                    auxnamespersource.append(sl[1])
                    data = np.array(alldatalist[isl],dtype=object)
                    dataperaux.append(data)
            # Indexed as uniquesources
            self.allauxnames.append( auxnamespersource )
            self.alldataperaux.append( dataperaux )
        
        # If mf6, also get the source kind for 
        # the given srcnames, just in case
        self.sourcestype = []
        if self.modelversion == 'mf6':
            for us in self.uniquesources:
                pkg = self._parent.flowmodel.get_package(pname.lower())
                self.sourcestype.append( pkg.package_type.upper() )

        # Save sources spec
        self.sources = sources

        # Define id
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER

        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = 'SRC'+str(self.__class__.COUNTER)


        # Add package and save instance        
        if self.__class__.COUNTER == 1:
            self.parent.add_package(self)
        self.__class__.INSTANCES.append( self )


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

        # Open file for writing
        f = open(self.INSTANCES[0].fn_path, "w")


        # Create format        
        fmts = []
        fmts.append("{:20s}")  # auxvarname
        fmts.append("{:.6f}")  # particles mass
        for nt in range(3):    # (nx,ny,nz)
            fmts.append("{:3d}")
        # Append solute id ?
        istart = 0
        iend   = 5
        if ( 
            ( self.INSTANCES[0]._parent.particlesmassoption == 2 ) or  # If solute id for particle groups
            ( self.INSTANCES[0]._parent.particlesmassoption == 2 ) ):  # Solute specific dispersion
            fmts.append("{:4d}")
            iend = 6
        fmt  = " ".join(fmts) + "\n"

        srcfmts = []
        srcfmts.append("{:20s}") # srcname
        srcfmt = " ".join(srcfmts) + "\n"

        # Write how many INSTANCES 
        # and loop over them
        f.write(f"{self.COUNTER}\n")

        for ins in self.INSTANCES:

            # Write id's
            f.write(f"{ins.stringid}\n")

            # Write format
            f.write(f"{ins.format.upper()}\n")

            if (
                ( ins.format.upper() == 'AUX' ) or 
                ( ins.format.upper() == 'AUXILIARY' ) ):

                # Inform about the number of source budgets
                f.write(f"{len(ins.uniquesources)}\n")

                # (...) and loop over them
                for isrc, src in enumerate(ins.uniquesources):

                    # Write src name
                    f.write(f"{src.upper()}\n")

                    # Write the number of aux vars
                    f.write(f"{len(ins.allauxnames[isrc])}\n")
                    
                    # Write aux vars
                    for iaux, auxdata in enumerate(ins.alldataperaux[isrc]):
                        f.write(fmt.format(*auxdata[istart:iend]))


        # And close
        f.close()


        # Done
        return
