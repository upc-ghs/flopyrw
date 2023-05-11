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









