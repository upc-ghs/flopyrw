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

    While writing input files with check=True, the existence of the 
    following packages will be verified: 
    
        - ModpathRWSim : the simulation class
        - ModpathRWDsp : dispersion package
        - ModpathRWOpts: random walk options 
        - ModpathRWBas : basic package, porosities
    '''

    def __init__(self, *args, **kwargs):

        # Assign defaults values, prioritizing user provided information
        try:
            kwargs['exe_name']
        except KeyError:
            kwargs['exe_name'] = 'mpathrw'
        try:
            kwargs['simfile_ext']
        except KeyError:
            kwargs['simfile_ext'] = 'mprw' # should be consistent with ModpathRWSim
        try:
            kwargs['modelname']
        except KeyError:
            kwargs['modelname'] = 'mprwsim' 


        # Call parent constructor
        super().__init__(*args,**kwargs)
        # update version
        self.version_types = {"modpathrw": "MODPATH-RW"}
        self.set_version('modpathrw')


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


    # Overload the write_input method
    def write_input(self, SelPackList=False, check=False):
        """
        Write the input.

        Parameters
        ----------
        SelPackList : False or list of packages

        Note: is practically the same than BaseModel but passing
        the check parameter until the write_name_file function
        """
        if check:
            # run check prior to writing input
            self.check(f=f"{self.name}.chk", verbose=self.verbose, level=1)

        # reset the model to free_format if parameter substitution was
        # performed on a model load
        if self.parameter_load and not self.free_format_input:
            if self.verbose:
                print(
                    "\nResetting free_format_input to True to "
                    "preserve the precision of the parameter data."
                )
            self.free_format_input = True

        if self.verbose:
            print("\nWriting packages:")

        if SelPackList == False:
            for p in self.packagelist:
                if self.verbose:
                    print("   Package: ", p.name[0])
                # prevent individual package checks from running after
                # model-level package check above
                # otherwise checks are run twice
                # or the model level check procedure would have to be split up
                # or each package would need a check argument,
                # or default for package level check would have to be False
                try:
                    p.write_file(check=False)
                except TypeError:
                    p.write_file()
        else:
            for pon in SelPackList:
                for i, p in enumerate(self.packagelist):
                    if pon in p.name:
                        if self.verbose:
                            print("   Package: ", p.name[0])
                        try:
                            p.write_file(check=False)
                        except TypeError:
                            p.write_file()
                            break
        if self.verbose:
            print(" ")

        # write name file
        self.write_name_file(check=check)
        # os.chdir(org_dir)



    # Overload the write_name_file method
    def write_name_file(self,check=False):
        """
        Write the name file
        
        check: bool
            Verifies consistency of the model before writing

        Returns
        -------
        None
        """


        with open(os.path.join(self.model_ws, self.mpnamefile), "w") as f:

            # Write the mpnam file
            f.write(f"{self.heading}\n")

            # check consistency 
            if check: 
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
                            f"{self.__class__.__name__}:"
                            f" A simulation package is necessary for modpath models,"
                            f" but None was given. Define a ModpathRWSim or Modpath7Sim package."
                        )

                # enforce bas
                baspkg = self.get_package('MPRWBAS')
                if baspkg is None:
                    baspkg = self.get_package('MPBAS')
                    if baspkg is None: 
                        raise Exception( 
                            f"{self.__class__.__name__}:"
                            f" The basic package is necessary for modpath models,"
                            f" but None was given. Define a ModpathRWBas or Modpath7Bas package."
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
            if check:
                if ( (dsppkg is None) and isrw ): 
                    raise Exception( 
                        f"{self.__class__.__name__}:" 
                        f" A dispersion package is required for a random walk model,"
                        f" but None was given. Define a ModpathRWDsp package."
                    )
            if dsppkg is not None:
                f.write(f"DSP        {dsppkg.file_name[0]}\n")
            # rwopts ( enforced if isrw )
            rwoptspkg = self.get_package('RWOPTS')
            if check:
                if ( (rwoptspkg is None) and isrw and check ): 
                    raise Exception( 
                        f"{self.__class__.__name__}:"
                        f" A random walk options package is required for a random walk model," 
                        f" but None was given. Define a ModpathRWOpts package." 
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
        return

