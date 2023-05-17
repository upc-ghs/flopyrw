'''
Configuration of MODPATH-RW 
'''

# python
import os

# flopy
from flopy.modpath import Modpath7


class ModpathRW( Modpath7 ):
    '''
    MODPATH-RW class 

    Extends flopy.modpath.Modpath7 and overloads method write_name_file

    A random walk simulation requires the DSP and RWOPTS packages
    and giving these will be enforced by the class before completing 
    the writing of input files.

    The basic package (porosities) is also required by modpath and this
    class will also enforce its specification if not given.

    A simulation class is also required by modpath and will be enforced 
    by the class.

    '''

    def __init__(
            self, *args,
            version='modpathrw',
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

        # Write the mpnam file
        f.write(f"{self.heading}\n")

        # enforce sim
        simpkg   = self.get_package('MPRWSIM')
        isrw     = False
        foundsim = False
        if simpkg is not None:
            foundsim = True
            if simpkg.simulationtype >= 5:
                isrw = True
        if not foundsim: 
            simpkg = self.get_package('MPSIM')
            if simpkg is None: 
                raise Exception( 
                    self.__class__.__name__ + ':' + 
                    ' A simulation package is necessary for modpath models, ' + 
                    'but None was given. Define a ModpathRWSim or Modpath7Sim package. ' 
                )

        # enforce bas
        baspkg = self.get_package('MPRWBAS')
        if baspkg is None:
            baspkg = self.get_package('MPBAS')
            if baspkg is None: 
                raise Exception( 
                    self.__class__.__name__ + ':' + 
                    ' The basic package is necessary for modpath models, ' + 
                    'but None was given. Define a ModpathRWBas or Modpath7Bas package. ' 
                )
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
        # dsp ( enforced if isrw )
        dsppkg = self.get_package('DSP')
        if ( (dsppkg is None) and isrw ): 
            raise Exception( 
                self.__class__.__name__ + ':' + 
                ' A dispersion package is required for a random walk model, ' + 
                'but None was given. Define a ModpathRWDsp package. ' 
            )
        if dsppkg is not None:
            f.write(f"DSP        {dsppkg.file_name[0]}\n")
        # rwopts ( enforced if isrw )
        rwoptspkg = self.get_package('RWOPTS')
        if ( (rwoptspkg is None) and isrw ): 
            raise Exception( 
                self.__class__.__name__ + ':' + 
                ' A random walk options package is required for a random walk model, ' + 
                'but None was given. Define a ModpathRWOpts package. ' 
            )
        if rwoptspkg is not None:
            f.write(f"RWOPTS     {rwoptspkg.file_name[0]}\n")
        # spc
        spcpkg = self.get_package('SPC')
        if spcpkg is not None:
            f.write(f"SPC        {spcpkg.file_name[0]}\n")
        # gpkde
        gpkdepkg = self.get_package('GPKDE')
        if gpkdepkg is not None:
            f.write(f"GPKDE      {gpkdepkg.file_name[0]}\n")
        # ic
        icpkg = self.get_package('IC')
        if icpkg is not None:
            f.write(f"IC         {icpkg.file_name[0]}\n")
        # imp
        imppkg = self.get_package('IMP')
        if imppkg is not None:
            f.write(f"IMP        {imppkg.file_name[0]}\n")
        # src
        srcpkg = self.get_package('SRC')
        if srcpkg is not None:
            f.write(f"SRC        {srcpkg.file_name[0]}\n")
        # obs
        obspkg = self.get_package('OBS')
        if obspkg is not None:
            f.write(f"OBS        {obspkg.file_name[0]}\n")


        # Done
        f.close()


        return
