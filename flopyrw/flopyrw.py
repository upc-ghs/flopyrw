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


# Utils
def count_instances(method):
    '''
    Util for counting class instances. 
    Rely on the class having a COUNTER variable
    '''
    def wrapper(self, *args, **kw):
        self.__class__.COUNTER += 1
        return method(self, *args, **kw)
    return wrapper


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



class ModpathRWDsp( Package ):
    """
    MODPATH-RW Dispersion Package Class.
    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    """


    # Class properties
    COUNTER    = 0
    UNITNUMBER = 0
    INSTANCES  = []


    @count_instances
    def __init__(
        self,
        model,
        modelkind         = None , # linear, nonlinear
        modelkindid       = None , # 1:linear, 2:nonlinear
        alphal            = 0.1  , # xx
        alphath           = 0.01 , # yy
        alphatv           = 0.01 , # zz
        dmeff             = 0.0  , # corrected by tortuosity
        dmaqueous         = 0.0  ,
        betal             = 1    ,
        betath            = 0.5  , 
        betatv            = 0.5  , 
        delta             = 5    , 
        dgrain            = 1    ,
        id                = None , 
        stringid          = None ,
        extension         = 'dsp',
    ):


        # Define UNITNUMBER if the first instance is created
        if self.__class__.COUNTER == 1:
            unitnumber = model.next_unit()
            super().__init__(model, extension, "DSP", unitnumber)
            self.__class__.UNITNUMBER = self.unit_number[0]
        else:
            # Needed for Util3d
            self._parent = self.INSTANCES[0]._parent
        

        # What about this ?
        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])
        self.model = model
        self.shape3d = shape3d


        # Select modelkind
        if modelkind is not None: 
            if modelkind not in ('linear', 'nonlinear'):
                raise Exception( 'flopyrwpt.py: modelkind ' + modelkind + ' is not allowed: linear or nonlinear' )
            if modelkind == 'linear':
                self.modelkind   = modelkind
                self.modelkindid = 1
            elif modelkind == 'nonlinear':
                self.modelkind   = modelkind
                self.modelkindid = 2
        else:
            self.modelkind   = 'linear'
            self.modelkindid = 1


        # Read dispersion model parameters 
        if self.modelkindid == 1:
            # linear
            self.alphal = Util3d(
                model,
                shape3d,
                np.float32,
                alphal,
                name="ALPHAL",
                locat=self.__class__.UNITNUMBER,
            )
            self.alphath = Util3d(
                model,
                shape3d,
                np.float32,
                alphath,
                name="ALPHATH",
                locat=self.__class__.UNITNUMBER,
            )
            self.alphatv = Util3d(
                model,
                shape3d,
                np.float32,
                alphatv,
                name="ALPHATV",
                locat=self.__class__.UNITNUMBER,
            )
            # Consider some checks to avoid writing unnecessary zeroes
            self.dmeff = Util3d(
                model,
                shape3d,
                np.float32,
                dmeff,
                name="DMEFF",
                locat=self.__class__.UNITNUMBER,
            )
        elif self.modelkindid == 2:
            # nonlinear
            self.dmaqueous = dmaqueous
            # Consider some checks to avoid writing unnecessary zeroes
            self.dmeff = Util3d(
                model,
                shape3d,
                np.float32,
                dmeff,
                name="DMEFF",
                locat=self.__class__.UNITNUMBER,
            )
            self.betal = Util3d(
                model,
                shape3d,
                np.float32,
                betal,
                name="BETAL",
                locat=self.__class__.UNITNUMBER,
            )
            self.betath= Util3d(
                model,
                shape3d,
                np.float32,
                betath,
                name="BETATH",
                locat=self.__class__.UNITNUMBER,
            )
            self.betatv= Util3d(
                model,
                shape3d,
                np.float32,
                betatv,
                name="BETATV",
                locat=self.__class__.UNITNUMBER,
            )
            self.delta = Util3d(
                model,
                shape3d,
                np.float32,
                delta,
                name="DELTA",
                locat=self.__class__.UNITNUMBER,
            )
            self.dgrain = Util3d(
                model,
                shape3d,
                np.float32,
                dgrain,
                name="DGRAIN",
                locat=self.__class__.UNITNUMBER,
            )



        # Define id
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER

        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = 'DSP'+str(self.__class__.COUNTER)


        # Add package and save instance        
        if self.__class__.COUNTER == 1:
            self.parent.add_package(self)
        self.__class__.INSTANCES.append( self )


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
            
            # Write id's
            f.write(f"{ins.id}\n")
            f.write(f"{ins.stringid}\n")

            # Write modelkindid
            f.write(f"{ins.modelkindid}\n")
           
            # Write dispersion model parameters
            if ins.modelkindid == 1:
                f.write(ins.dmeff.get_file_entry())
                f.write(ins.alphal.get_file_entry())
                f.write(ins.alphath.get_file_entry())
                f.write(ins.alphatv.get_file_entry())

            if ins.modelkindid == 2:
                f.write(f"{ins.dmaqueous:.16f}\n")
                f.write(ins.dmeff.get_file_entry())
                f.write(ins.betal.get_file_entry())
                f.write(ins.betath.get_file_entry())
                f.write(ins.betatv.get_file_entry())
                f.write(ins.delta.get_file_entry())
                f.write(ins.dgrain.get_file_entry())
       

        # Done
        f.close()

        return




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



class ModpathRWObs( Package ):
    """
    MODPATH-RW Observation Package Class.
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
        kind           = 1,
        celloption     = 1,
        cells          = None,
        timeoption     = 1,
        structured     = True, 
        basefilename   = 'mpathrwobs_',
        filename       = None,
        id             = None,
        stringid       = None,
        extension      = 'obs',
    ):
       

        # Define UNITNUMBER if the first instance is created
        if self.__class__.COUNTER == 1:
            unitnumber = model.next_unit()
            super().__init__(model, extension, "OBS", unitnumber)
            self.__class__.UNITNUMBER = self.unit_number[0]
        else:
            # Needed for Util3d
            self._parent = self.INSTANCES[0]._parent


        # Need health checks
        self.kind = kind
        self.timeoption = timeoption
        self.structured = structured

        # celloption 
        if ( celloption not in [1,2] ):
            raise ValueError('Invalid celloption ',
                    celloption, '. Allowed values are 1 (list of cellids) or 2 (modelgrid like array)')
        self.celloption = celloption

        # Write as a list of cellids
        if self.celloption == 1:
            if cells is None:
                self.cells = []
            else:
                if not isinstance( cells, (list,np.ndarray,tuple) ):
                    raise Exception( 'Cells parameter should be a list or numpy array, is ', type( cells ))
                if ( isinstance( cells, tuple ) and len(cells)==1 ):
                    cells = cells[0]
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
                    locat=self.__class__.UNITNUMBER,
                )


        # Define obs id
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER

        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = 'OBS'+str(self.__class__.COUNTER)

        # Filename for this observation
        self.filename = basefilename+str(self.id)+'.'+extension

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

            # Write obs id
            f.write(f"{ins.id}\n")

            # Write obs filename
            f.write(f"{ins.filename}\n")

            # Write the obs kind
            f.write(f"{ins.kind}\n")

            # Write the celloption param
            f.write(f"{ins.celloption}\n")

            if ins.celloption == 1:
                # Should write a cell number
                if len( ins.cells ) == 0: 
                    raise Exception('flopyrw:ModpathRWObs: cells array is empty. Specify a list of cells for the observation id ', self.id )
                else:
                    f.write(f"{len(ins.cells)}\n")
                    fmts = []
                    if ins.structured:
                        f.write(f"1\n") # To indicate structured
                        fmts.append("{:9d}") # lay
                        fmts.append("{:9d}") # row
                        fmts.append("{:9d}") # col
                    else:
                        f.write(f"2\n") # To indicate cell ids
                        fmts.append("{:9d}") # cellid
                    fmt = " " + " ".join(fmts) + "\n"
                    for oc in ins.cells:
                        woc = np.array(oc).astype(np.int32)+1 # Correct the zero-based indexes
                        if ins.structured:
                            f.write(fmt.format(*woc))
                        else:
                            f.write(fmt.format(woc))

            elif ins.celloption == 2:
                # Should write an array
                # with the distribution of the observation
                f.write(ins.cells.get_file_entry())


            # Write timeoption params
            f.write(f"{ins.timeoption}\n")
            if ins.timeoption == 1:
                pass
            elif ins.timeoption == 2:
                pass


        # Done
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


class ModpathRWBc( Package ):
    """
    MODPATH-RW Boundary Conditions Package Class.
    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    icbound : int, np.ndarray
        Determine rebound boundary conditions for RW particles.
    """


    def __init__(
        self,
        model,
        icbound    = None,
        flux       = None,
        prescribed = None, 
        extension  = 'bc',
    ):
       
        # Define UNITNUMBER if the first instance is created
        unitnumber = model.next_unit()
        super().__init__(model, extension, "BC", unitnumber)

        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])
        self.model = model
        self.shape3d = shape3d

        # icbound
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
            # Particles will not rebound
            self.icbound= Util3d(
                model,
                shape3d,
                np.int32,
                0,
                name="ICBOUND",
                locat=self.unit_number[0],
            )

        # flux
        if flux is not None:
            if isinstance(flux,list):
                for pic in flux:
                    if not isinstance(pic, ModpathRWFlux):
                        raise TypeError('flopyrw:ModpathRWBc: Object in list is type ', type(pic), '. Expected ModpathRWFlux.')
                # Survived so continue
                self.npfs = len(flux)
                self.flux = flux
            elif isinstance( flux, ModpathRWFlux ):
                self.npfs = 1
                self.flux = [flux]
            else:
                raise TypeError('flopyrw:ModpathRWBc: Flux argument should be of type list or ModpathRWFlux. ', type(flux), ' given.')
        else:
            self.npfs = 0
            self.flux = None


        # prescribed concentrations
        if prescribed is not None:
            if isinstance(prescribed,list):
                for pic in prescribed:
                    if not isinstance(pic, ModpathRWPc):
                        raise TypeError('flopyrw:ModpathRWBc: Object in list is type ', type(pic), '. Expected ModpathRWPc.')
                # Survived so continue
                self.npcs = len(prescribed)
                self.prescribed = prescribed
            elif isinstance( prescribed, ModpathRWPc ):
                self.npcs = 1
                self.prescribed = [prescribed]
            else:
                raise TypeError('flopyrw:ModpathRWBc: Prescribed argument should be of type list or ModpathRWPc. ', type(prescribed), ' given.')
        else:
            self.npcs = 0
            self.prescribed = None


        # Add package 
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


        # icbound 
        f.write(self.icbound.get_file_entry())

        # REQUIRES particlesmassoption ?
        # REQUIRES solutesoption ?


        # flux
        if( ( self.flux is None ) or ( self.npfs==0 ) ):
            # no flux
            f.write(f"0\n")
        else:
            # write fluxes
            f.write(f"{self.npfs}\n")

            for ins in self.flux:
                ins.write(
                  f=f,
                  #particlesmassoption=self.particlesmassoption,
                  #solutesoption=self.solutesoption,
                )


        # prescribed
        if( ( self.prescribed is None ) or ( self.npcs==0 ) ):
            # no prescribed
            f.write(f"0\n")
        else:
            # write prescribed
            f.write(f"{self.npcs}\n")

            for ins in self.prescribed:
                ins.write(
                  f=f,
                  #particlesmassoption=self.particlesmassoption,
                  #solutesoption=self.solutesoption,
                )


        # And close
        f.close()

        # Done
        return







class ModpathRWPc( Package ):
    """
    MODPATH-RW Prescribed Concentration Class.
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
        particlesmass = 0.1,
        concentration = 1.0,
        kind          = 1,
        celloption    = 1,
        structured    = True,
        cells         = None,
        id            = None,
        stringid      = None,
    ):
       

        ## Define UNITNUMBER if the first instance is created
        #if self.__class__.COUNTER == 1:
        #    unitnumber = model.next_unit()
        #    super().__init__(model, extension, "OBS", unitnumber)
        #    self.__class__.UNITNUMBER = self.unit_number[0]
        #else:
        #    # Needed for Util3d
        #    self._parent = self.INSTANCES[0]._parent


        self.particlesmass = particlesmass
        
        self.concentration = concentration




        # Need health checks
        self.kind = kind
        #self.timeoption = timeoption
        self.structured = structured

        # celloption 
        if ( celloption not in [1,2] ):
            raise ValueError('Invalid celloption ',
                    celloption, '. Allowed values are 1 (list of cellids) or 2 (modelgrid like array)')
        self.celloption = celloption

        # Write as a list of cellids
        if self.celloption == 1:
            if cells is None:
                self.cells = []
            else:
                if not isinstance( cells, (list,np.ndarray,tuple) ):
                    raise Exception( 'Cells parameter should be a list or numpy array, is ', type( cells ))
                if ( isinstance( cells, tuple ) and len(cells)==1 ):
                    cells = cells[0]
                # Maybe some sanity check about data structure or the same 
                # used for partlocs
                self.cells = cells

        # Write obs cells as 3D array
        if self.celloption == 2: 
            
            # NOT YET
            pass

            ## This was already done right?
            #shape = model.shape
            #if len(shape) == 3:
            #    shape3d = shape
            #elif len(shape) == 2:
            #    shape3d = (shape[0], 1, shape[1])
            #else:
            #    shape3d = (1, 1, shape[0])
            #self.model   = model
            #self.shape3d = shape3d

            #if cells is None:
            #    self.cells = []
            #else:
            #    if not isinstance( cells, (list,np.ndarray) ):
            #        raise Exception( 'Cells parameter should be a list or numpy array')
            #    self.cells = Util3d(
            #        model,
            #        shape3d,
            #        np.int32,
            #        cells,
            #        name="PCELLS",
            #        locat=self.__class__.UNITNUMBER,
            #    )


        # Define id
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER

        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = 'PC'+str(self.__class__.COUNTER)

        # Filename for this observation
        #self.filename = basefilename+str(self.id)+'.'+extension

        # Add package and save instance        
        #if self.__class__.COUNTER == 1:
        #    self.parent.add_package(self)
        self.__class__.INSTANCES.append( self )


        # Done
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
            raise Exception('ModpathRWPc: requires file pointer f. Is None.')

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


        # Needs some kind of soluteid


        # Needs the particlesmass
        f.write(f"{self.particlesmass:.10f}\n")


        # And the prescribed concentration
        f.write(f"{self.concentration:.10f}\n")



        # Write the celloption param
        f.write(f"{self.celloption}\n")

        if self.celloption == 1:
            # Should write a cell number
            if len( self.cells ) == 0: 
                raise Exception('flopyrw:ModpathRWPc: cells array is empty. Specify a list of cells for the observation id ', self.id )
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
                    if self.structured:
                        f.write(fmt.format(*woc))
                    else:
                        f.write(fmt.format(woc))

        elif self.celloption == 2:
            print( 'NOT IMPLEMENTED! ')
            pass
            ## Should write an array
            ## with the distribution of the observation
            #f.write(self.cells.get_file_entry())


        ## Write timeoption params
        #f.write(f"{self.timeoption}\n")
        #if self.timeoption == 1:
        #    pass
        #elif self.timeoption == 2:
        #    pass


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


































#########################################################
# DEPRECATION WARNING
#########################################################

class ModpathRWFlux( Package ): # Ssm ?
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
        super().__init__(model, extension, "FLUX", unitnumber)

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


# DEPRECATION WARNING
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
                if not isinstance( cells, (list,np.ndarray,tuple) ):
                    raise Exception( 'Cells parameter should be a list or numpy array, is ', type( cells ))
                if ( isinstance( cells, tuple ) and len(cells)==1 ):
                    cells = cells[0]
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
                    if self.structured:
                        f.write(fmt.format(*woc))
                    else:
                        f.write(fmt.format(woc))

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
        #if reconstruction:
        #    if reconstructionfilename is None:
        #        reconstructionfilename = f"{model.name}.gpkde"
        #    self.reconstructionfilename = reconstructionfilename
        #self.reconstruction = reconstruction

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
        super().__init__(model, extension, "FLUX", unitnumber)

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




# DEPRECATION WARNING ?
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
