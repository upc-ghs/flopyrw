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


class timeseriesOutputOption(Enum): 
    '''
    Enumerate formats for timeseries output option
    '''
    active = 0
    all    = 1
    skip   = 2


class endpointOutputOption(Enum): 
    '''
    Enumerate formats for endpoint output option 
    '''
    classic = 0
    mass    = 1
    skip    = 2


class particlesMassOption(Enum): 
    '''
    Enumerate formats for particles mass option 
    '''
    off     = 0
    mass    = 1
    massspc = 2


class speciesDispersionOption(Enum): 
    '''
    Enumerate formats for species dispersion option 
    '''
    unique   = 0
    specific = 1


class ModpathRWSim( mp7.Modpath7Sim ):
    '''
    MODPATH-RW Simulation File Package Class. 

    Extends from Modpath7Sim

    New Parameters
    --------------
    timeseriesoutputoption : int/ str
        Specify the behavior of timeseries writer. Allowed values are: 
         * 0 or 'active': Write timeseries records only for active particles.
         * 1 or 'all': Write timesereis records for all particles.
         * 2 or 'skip': Skip the timeseries writer.
    endpointoutputoption : int/ str
        Specify the behavior of endpoint writer. Allowed values are: 
         * 0 or 'classic': Write endpoint file with the classical format, modpath-v7. 
         * 1 or 'mass': Write endpoint file with the simplified format for mass particles. 
         * 2 or 'skip': Skip the endpoint writer.
    particlesmassoption: int/str
        Indicates whether a mass shall be read for classical particle groups
        specifications. Allowed values are: 
         * 0 or 'off': do not read a particle mass 
         * 1 or 'mass': read the particle mass for particle groups
         * 2 or 'massspc': read a mass and a solute identifier for all specified particle groups
    speciesdispersionoption: int/str
        Configures if particles are displaced with the same dispersion parameters or
        with species specific properties depending on solute id. Allowed values are:
         * 0 or 'unique': all particles displaced with the same dispersion properties
         * 1 or 'specific': particles displaced with specific dispersion properties
    '''

    
    @staticmethod
    def _ftype():
        return 'MPRWSIM'


    def __init__(
            self,
            *args,
            timeseriesoutputoption=0,
            endpointoutputoption=0,
            particlesmassoption=0,
            speciesdispersionoption=0,
            **kwargs
        ):

        # try to process an extension if given
        try:
            extension = kwargs['extension']
        except KeyError:
        # if not, extract from parent model
            extension = args[0].namefile_ext
        finally:
            extension = "mprw"

        # Call parent constructor
        # This function calls add_package().
        # It maybe worth considering to remove 
        # the package in case that any of exceptions 
        # below is triggered, to avoid manually removing
        # while testing.
        super().__init__(*args,**kwargs, extension=extension)
        
        # Not the most elegant solution but consistent
        ftype = self._ftype()
        self.name  = [ftype]
        self._name = [ftype]
        self._generate_heading()

        # Extract model
        model = args[0]

        # Override interpretation of particlegroups, allows None
        if 'particlegroups' not in kwargs.keys():
            particlegroups = None
            self.particlegroups = particlegroups

        # timeseriesoutputoption 
        if ( not isinstance( timeseriesoutputoption, (int,str) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for timeseriesoutputoption."
                f" It can be specified as int or str, but {str(type(timeseriesoutputoption))}"
                f" was given." 
            )
        if ( isinstance( timeseriesoutputoption, int ) ):
            if (timeseriesoutputoption not in [0,1,2]):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid timeseriesoutputoption {str(timeseriesoutputoption)}"
                    f" Allowed values are 0 (write ts for active particles),"
                    f" 1 (write ts for all particles) or 2 (skip ts writer)."
                )
            self.timeseriesoutputoption = timeseriesoutputoption
        elif ( isinstance( timeseriesoutputoption, str ) ):
            try:
                self.timeseriesoutputoption = timeseriesOutputOption[timeseriesoutputoption.lower()].value
            except:
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid timeseriesoutputoption {str(timeseriesoutputoption)}."
                    f" The allowed values are active, all or skip."
                )

        # endpointoutputoption
        if ( not isinstance( endpointoutputoption, (int,str) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for endpointoutputoption."
                f" It can be specified as int or str, but {str(type(endpointoutputoption))}"
                f" was given." 
            )
        if ( isinstance( endpointoutputoption, int ) ):
            if (endpointoutputoption not in [0,1,2]):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid endpointoutputoption {str(endpointoutputoption)}"
                    f" Allowed values are 0 (classic endpoint output),"
                    f" 1 (mass particles endpoint output) or "
                    f" 2 (skip writing of endpoint file)."
                )
            self.endpointoutputoption = endpointoutputoption
        elif ( isinstance( endpointoutputoption, str ) ):
            try:
                self.endpointoutputoption = endpointOutputOption[endpointoutputoption.lower()].value
            except:
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid endpointoutputoption {str(endpointoutputoption)}."
                    f" The allowed values are classic or mass."
                )

        # particlesmassoption
        if ( not isinstance( particlesmassoption, (int,str) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for particlesmassoption."
                f" It can be specified as int or str, but {str(type(particlesmassoption))}"
                f" was given." 
            )
        if ( isinstance( particlesmassoption, int ) ):
            if (particlesmassoption not in [0,1,2]):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid particlesmassoption {str(particlesmassoption)}."
                    f" Allowed values are 0 (no particle mass), "
                    f" 1 (read mass) or 2 (read mass and a solute id)."
                )
            self.particlesmassoption = particlesmassoption
        elif ( isinstance( particlesmassoption, str ) ):
            try:
                self.particlesmassoption = particlesMassOption[particlesmassoption.lower()].value
            except:
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid particlesmassoption {str(particlesmassoption)}."
                    f" The allowed values are off, mass or massspc."
                )

        # speciesdispersionoption
        if ( not isinstance( speciesdispersionoption, (int,str) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for speciesdispersionoption."
                f" It can be specified as int or str, but {str(type(speciesdispersionoption))}"
                f" was given." 
            )
        if ( isinstance( speciesdispersionoption, int ) ):
            if (speciesdispersionoption not in [0,1]):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid speciesdispersionoption {str(speciesdispersionoption)}"
                    f" Allowed values are 0 (same dispersion parameters) or"
                    f" 1 (species specific dispersion)."
                )
            self.speciesdispersionoption = speciesdispersionoption
        elif ( isinstance( speciesdispersionoption, str ) ):
            try:
                self.speciesdispersionoption = speciesDispersionOption[speciesdispersionoption.lower()].value
            except:
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid speciesdispersionoption {str(speciesdispersionoption)}."
                    f" The allowed values are unique or specific."
                )

        # If simulation is RW
        if self.simulationtype > 4: 
            # Assign some properties to parent obj
            # Needed by: SPC, IC, SRC
            self._parent.particlesmassoption = self.particlesmassoption
            self._parent.speciesdispersionoption = self.speciesdispersionoption

        
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

        Note: mostly the same as classic modpath but with management of rwpt options
        """
    
        with open(self.fn_path, "w") as f:

            # item 0
            f.write(f"{self.heading}\n")
            # item 1
            f.write(f"{self.mp_name_file}\n")
            # item 2
            f.write(f"{self.listingfilename}\n")
            # item 3
            f.write(
                "{} {} {} {} {} {} {} {} {} {}\n".format(
                    self.simulationtype,
                    self.trackingdirection,
                    self.weaksinkoption,
                    self.weaksourceoption,
                    self.budgetoutputoption,
                    self.tracemode,
                    self.timeseriesoutputoption,
                    self.endpointoutputoption,
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
 

        # Done
        return
