'''
Configuration of MODPATH-RW 
'''

# python
import os
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen

# flopy
from flopy.modpath import Modpath7
from typing import List, Optional, Tuple, Union
from flopy.mbase import resolve_exe
from flopy.utils import flopy_io


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
        try:
            kwargs['flowmodel']
            # If flowmodel was given and not model_ws, 
            # define by default with the model_ws from 
            # the flowmodel
            try:
                kwargs['model_ws']
            except KeyError:
                kwargs['model_ws'] = kwargs['flowmodel'].model_ws
        except KeyError:
            raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" The flowmodel argument should be given for defining a particle tracking model."
                )

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


    def run_model(
        self,
        silent=False,
        pause=False,
        report=False,
        normal_msg="normal termination",
        nprocesses=None,
    ) -> Tuple[bool, List[str]]:
        """
        This method will run the model using subprocess.Popen.
        It is based on the original run_model function from BaseModel 
        in flopy. Some additional logic is included in order to handle
        specific use cases from MODPATH-RW.

        Parameters
        ----------
        silent : boolean
            Echo run information to screen (default is True).
        pause : boolean, optional
            Pause upon completion (default is False).
        report : boolean, optional
            Save stdout lines to a list (buff) which is returned
            by the method . (default is False).
        normal_msg : str
            Normal termination message used to determine if the
            run terminated normally. (default is 'normal termination')
        nprocesses : int
            Defines the number of threads to be used in parallel 
            simulations. It should be greater than 0.
    
        Returns
        -------
        success : boolean
        buff : list of lines of stdout
    
        """
   

        # Interpret nprocesses
        cargs = None
        if nprocesses is not None:
            if not isinstance(nprocesses,int):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid nprocesses parameter of type {str(type(nprocesses))}."
                    f" It should be of type integer."
                )
            if nprocesses <= 0:
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid nprocesses, it should be greater than zero."
                )

            # If valid nprocesses, then format a cargs
            cargs = ["-np", f"{str(nprocesses)}"]
       
        return run_model(
            self.exe_name,
            self.namefile,
            model_ws=self.model_ws,
            silent=silent,
            pause=pause,
            report=report,
            normal_msg=normal_msg,
            cargs=cargs,
        )


def run_model(
    exe_name: Union[str, os.PathLike],
    namefile: Optional[str],
    model_ws: Union[str, os.PathLike] = os.curdir,
    silent=False,
    pause=False,
    report=False,
    processors=None,
    normal_msg="normal termination",
    use_async=False,
    cargs=None,
) -> Tuple[bool, List[str]]:
    """
    This function is mostly the same as the original from flopy, the only 
    difference being that the command arguments (cargs) are passed before the namefile,
    which is the structure in MODPATH-RW.
    
    Original documentation
    ----------------------

    Run the model using subprocess.Popen, optionally collecting stdout and printing
    timestamped progress. Model workspace, namefile, executable to use, and several
    other options may be configured, and additional command line arguments may also
    be provided.

    Parameters
    ----------
    exe_name : str or PathLike
        Executable name or path. If the executable name is provided,
        the executable must be on the system path. Alternatively, a
        full path to the executable may be provided.
    namefile : str, optional
        Name of the name file of model to run. The name may be None
        to run models that don't require a control file (name file)
    model_ws : str or PathLike, optional, default '.'
        Path to the parent directory of the namefile. (default is the
        current working directory '.')
    silent : boolean, default True
        Whether to suppress model output. (Default is True)
    pause : boolean, optional, default False
        Pause and wait for keystroke upon completion. (Default is False)
    report : boolean, optional, default False
        Save stdout lines to a list (buff) returned by the method. (Default is False)
    processors: int
        Number of processors. Parallel simulations are only supported for
        MODFLOW 6 simulations. (default is None)
    normal_msg : str or list
        Termination message used to determine if the model terminated normally.
        More than one message can be provided using a list.
        (Default is 'normal termination')
    use_async : boolean
        Asynchronously read model stdout and report with timestamps. Good for
        models taking a long time to run, not good for models that run quickly.
    cargs : str or list, optional, default None
        Additional command line arguments to pass to the executable.
        (Default is None)
    Returns
    -------
    success : boolean
    buff : list of lines of stdout (empty if report is False)

    """
    success = False
    buff = []

    # convert normal_msg to a list of lower case str for comparison
    if isinstance(normal_msg, str):
        normal_msg = [normal_msg]
    for idx, s in enumerate(normal_msg):
        normal_msg[idx] = s.lower()

    # make sure executable exists
    if exe_name is None:
        raise ValueError(f"An executable name or path must be provided")
    exe_path = resolve_exe(exe_name)
    if not silent:
        print(
            f"FloPy is using the following executable to run the model: {flopy_io.relpath_safe(exe_path, model_ws)}"
        )

    # make sure namefile exists
    if namefile is not None and not os.path.isfile(
        os.path.join(model_ws, namefile)
    ):
        raise FileNotFoundError(
            f"The namefile for this model does not exist: {namefile}"
        )

    # simple little function for the thread to target
    def q_output(output, q):
        for line in iter(output.readline, b""):
            q.put(line)
            # time.sleep(1)
            # output.close()

    # create a list of arguments to pass to Popen
    if processors is not None:
        if "mf6" not in exe_path:
            raise ValueError("processors kwarg only supported for MODFLOW 6")
        mpiexec_path = resolve_exe("mpiexec")
        if not silent:
            print(
                f"FloPy is using {mpiexec_path} "
                + f"to run {exe_path} "
                + f"on {processors} processors."
            )
        argv = [mpiexec_path, "-np", f"{processors}", exe_path, "-p"]
    else:
        argv = [exe_path]

    # This is the difference with the original function in flopy.mbase 

    # add additional arguments to Popen arguments
    if cargs is not None:
        if isinstance(cargs, str):
            cargs = [cargs]
        for t in cargs:
            argv.append(t)

    if namefile is not None:
        argv.append(Path(namefile).name)

    # END This is the difference with the original function in flopy.mbase 

    # run the model with Popen
    proc = Popen(argv, stdout=PIPE, stderr=STDOUT, cwd=model_ws)

    if not use_async:
        while True:
            line = proc.stdout.readline().decode("utf-8")
            if line == "" and proc.poll() is not None:
                break
            if line:
                for msg in normal_msg:
                    if msg in line.lower():
                        success = True
                        break
                line = line.rstrip("\r\n")
                if not silent:
                    print(line)
                if report:
                    buff.append(line)
            else:
                break
        return success, buff

    # some tricks for the async stdout reading
    q = Queue.Queue()
    thread = threading.Thread(target=q_output, args=(proc.stdout, q))
    thread.daemon = True
    thread.start()

    failed_words = ["fail", "error"]
    last = datetime.now()
    lastsec = 0.0
    while True:
        try:
            line = q.get_nowait()
        except Queue.Empty:
            pass
        else:
            if line == "":
                break
            line = line.decode().lower().strip()
            if line != "":
                now = datetime.now()
                dt = now - last
                tsecs = dt.total_seconds() - lastsec
                line = f"(elapsed:{tsecs})-->{line}"
                lastsec = tsecs + lastsec
                buff.append(line)
                if not silent:
                    print(line)
                for fword in failed_words:
                    if fword in line:
                        success = False
                        break
        if proc.poll() is not None:
            break
    proc.wait()
    thread.join(timeout=1)
    buff.extend(proc.stdout.readlines())
    proc.stdout.close()

    for line in buff:
        for msg in normal_msg:
            if msg in line.lower():
                print("success")
                success = True
                break

    if pause:
        input("Press Enter to continue...")
    return success, buff

