'''
Configuration of MODPATH-RW particle groups

Update particle groups classes from modpath7.
Extend interfaces to understand mass and solute
'''

# flopy
from flopy.modpath.mp7particlegroup import (
    ParticleGroup,
    ParticleGroupLRCTemplate,
    ParticleGroupNodeTemplate,
)
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
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Given mass for the particle group is zero. '
            )
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


class ParticleGroupLRCTemplate( mp7ParticleGroupLRCTemplate ): 

    def __init__( self, *args, mass=1.0, solute=0, **kwargs ):
        # Call parent constructor
        super().__init__(*args,**kwargs) 
        
        if mass != 0:
            self.mass = mass
        else:
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Given mass for the particle group is zero. '
            )
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


class ParticleGroupNodeTemplate( mp7ParticleGroupNodeTemplate ): 

    def __init__( self, *args, mass=1.0, solute=0, **kwargs ):
        # Call parent constructor
        super().__init__(*args,**kwargs) 
        
        if mass != 0:
            self.mass = mass
        else:
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Given mass for the particle group is zero. '
            )
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
