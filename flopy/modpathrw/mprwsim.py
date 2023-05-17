'''
Configuration of MODPATH-RW simulation
'''

# python
import numpy as np
from enum import Enum

# flopy
import flopy.modpath.mp7sim as mp7
from flopy.utils import Util2d


# Updates enumeration of available simulations
class simType(Enum):
    '''
    Enumeration of different simulation types
    
    With respect to modpath-v7, add:
        rwtimeseries
        rwcombined
        rwendpoint
    '''

    endpoint     = 1
    pathline     = 2
    timeseries   = 3
    combined     = 4
    rwtimeseries = 5
    rwcombined   = 6
    rwendpoint   = 7
mp7.simType = simType


class ModpathRWSim( mp7.Modpath7Sim ):
    '''
    MODPATH-RW Simulation File Package Class. 

    Extends from Modpath7Sim

    New Parameters
    --------------
    timeseriesoutputoption : int
        Specify behavior of timeseries writer. Allowed values are: 
            * 0: Write timeseries records only for active particles
            * 1: Write timesereis records for all particles
            * 2: Skip the timeseries writer
    particlesmassoption: int
        Indicates whether a mass shall be read for classical particle groups
        specifications. Allowed values are: 
            * 0: do not read a particle mass 
            * 1: read the particle mass for particle groups
            * 2: read a mass and a solute identifier for all specified particle groups
    speciesdispersionoption: int
        Configures if particles are displaced with the same dispersion parameters or
        with species specific properties depending on solute id. Allowed values are:
            * 0: all particles displaced with the same dispersion properties
            * 1: particles displaced with specific dispersion properties
    '''

    def __init__(
            self,
            *args,
            timeseriesoutputoption=0,
            particlesmassoption=0,
            speciesdispersionoption=0,
            **kwargs
        ):

        # Call parent constructor
        super().__init__(*args,**kwargs, extension='mprw')
        
        # Not the most elegant solution but consistent
        self.name  = ['MPRWSIM']
        self._name = ['MPRWSIM']
        self._generate_heading()

        # Extract model
        model = args[0]

        # Override interpretation of particlegroups, allows None
        if 'particlegroups' not in kwargs.keys():
            particlegroups = None
            self.particlegroups = particlegroups

        # New options
        if (timeseriesoutputoption not in [0,1,2]):
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid timeseriesoutputoption ' +str(timeseriesoutputoption)+
                '. Allowed values are 0 (write ts for active particles), ' +
                '1 (write ts for all particles) or 2 (skip ts writer).'
            )
        self.timeseriesoutputoption = timeseriesoutputoption
        if (particlesmassoption not in [0,1,2]):
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid particlesmassoption ' +str(particlesmassoption)+
                '. Allowed values are 0 (no particle mass), ' +
                '1 (read mass) or 2 (read mass and a solute id).'
            )
        self.particlesmassoption = particlesmassoption
        if (speciesdispersionoption not in [0,1]):
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid speciesdispersionoption ' + str(speciesdispersionoption)+
                '. Allowed values are 0 (same dispersion parameters) or ' +
                '1 (species specific dispersion).'
            )
        self.speciesdispersionoption = speciesdispersionoption


        # If simulation is RW
        if self.simulationtype > 4: 
            # Assign some properties to parent obj
            # Needed by: SPC, IC, SRC
            self._parent.particlesmassoption = particlesmassoption
            self._parent.speciesdispersionoption = speciesdispersionoption


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
                self.speciesdispersionoption,
            )
        )
        # item 4
        f.write(f"{self.endpointfilename}\n")
        # item 5
        if (
            self.simulationtype == 2 or
            self.simulationtype == 4 or
            self.simulationtype == 6
        ):
            f.write(f"{self.pathlinefilename}\n")
        # item 6
        if (
            self.simulationtype == 3 or
            self.simulationtype == 4 or
            self.simulationtype == 5 or
            self.simulationtype == 6
        ):
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
            self.simulationtype == 6
        ):
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
