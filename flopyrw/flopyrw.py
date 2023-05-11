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

