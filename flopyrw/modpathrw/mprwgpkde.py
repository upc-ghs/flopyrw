'''
Configuration of MODPATH-RW gpkde reconstruction
'''

# python
import numpy as np

# flopy
from flopy.pakbase import Package


class ModpathRWGpkde( Package ):
    """
    MODPATH-RW GPKDE Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class: ModpathRW) to which
        this package will be added.
    outputfilename: str
        Filename that will be used for writing reconstruction output.
    outputcolformat : int
        Determines the columns structure of the output file. Allowed values are: 
          * 0: bin ids and density data
          * 1: bin ids, cell coordinates and density data
          * 2: cell coordinates and density data
    outputfileformat : int
        Determines the format of the output file. Allowed values are: 
          * 0: text plain file
          * 1: binary unformatted file
    domainorigin : list [ox,oy,oz]
        List with domain origin. If not given will be set as default to 
        [ xoffset, yoffset, min(botm) ]
    domainsize : list [sx,sy,sz]
        List with domain sizes. These are used for computing reconstruction 
        grid dimensions using given binsize.
    binsize : list [sx,sy,sz]
        List with cell sizes employed for reconstruction process. If the size on 
        a given dimension is zero, then the histogram count will not consider said 
        dimension.
    gridallocformat : int
        Determines the allocation in memory of the reconstruction grid. Allowed values are: 
          * 0: allocate following domain size 
          * 1: allocate according to the particle distribution
    gridborderfraction : float
        While allocating according to the particle locations, border fraction determines 
        an additional allocated space surrounding the particle cloud, extends the 
        bounding box allowing the boundaries to receive information from the interior
        histogram cells. For a given dimension, the distance considered as buffer
        is equal to gridborderfraction*extentdimension and is distributed as half
        for each side. Should be between [0,1].
    noptloops : int
        The number of bandwidth optimization loops. It should be a positive integer.
    skiperror : bool
        Flag to indicate whether the error checks for the optimization of bandwidth
        selection should be skipped or not. If skipped, optimization is performed 
        until noptloops.
    convergence : float
        The limit relative change of the variables monitored during bandwidth optimization 
        to determine the convergence of reconstruction. It should be positive.
    kerneldatabase : bool
        Flag to indicate whether reconstruction should be performed with a database of 
        kernels (True) or by computing kernels in real time (False). In case a kernel 
        database is allocated, then the class will employ the parameters 
        minhd: deltahd: maxhd to determine a range of discrete kernel sizes for
        allocation. This range is relative to the bin size, meaning, 
        parameters *hd indicate the ratio smoothing(h)/binsize(d) for each model dimension.
    isotropickernels : bool
        Flag to indicate that kernels are considered as isotropic. This consideration 
        is done in terms of the kernel smoothing (h), meaning that if the bin sizes 
        are non-isotropic, then the program will force the same kernel smoothing taking 
        into account the different cell sizes.
    boundkernelsize : int
        Determines the protocol for bounding the kernel size. Allowed values are: 
          * 0: bound size by domain constraints. 
          * 1: bound using minhd and maxhd
          * 2: unbounded
    minhd : float
        Minimum ratio smoothing(h)/binsize(l). It is written to the package file 
        if kerneldatabase=True or boundkernelsize = 1.
    deltahd : float
        Step of ration smoothing(h)/binsize(l) used for defining a set of discrete
        kernel sizes for the kernel database. Only written to the package file 
        if kerneldatabase=True.
    maxhd : float
        Maximum ration smoothing(h)/binsize(l). It is written to the package file if
        kerneldatabase=True or boundkernelsize = 1.
    initialsmoothingformat : int
        Determines the protocol for selection of the initial kernel size from where to
        begin the optimization for bandwidth selection. Allowed values are: 
          * 0: selects initial bandwidth from the expression of Silverman (1986) 
               for Gaussian distributions. It is based on the standard deviation of
               the particle coordinates, and the number of particles. 
          * 1: initial bandwidth as a factor multiplying the characteristic cell size.
               It will use the parameter binsizefactor to perform the scaling.
    binsizefactor : float
        Scaling factor employed for determining the initial kernel size from where to start
        the bandwidth optimization. Amplifies the characteristic cell size.
    asconcentration: bool
        Flag to indicate whether the density reconstruction should be returned as 
        resident concentration (dissolved). Because the reconstruction grid is different
        than the flowmodel grid, said transformation can be easily achieved only in case
        both porosities and delaying factors are spatially uniform. Otherwise, density is 
        returned as total mass concentration, including both sorbed and aqueous phases.
    effectiveweightformat : int
        Determines the protocol for handling the reconstruction of particle distributions 
        with different weights. Allowed values are: 
          * 0: computes a unique effective number of particles and determines a 
               domain-effective particle weight (Kish, 1965,1992), which is used to
               transform the mass histogram into an equivalent particle distribution 
               used for the bandwidth selection.
          * 1: similar to the previous alternative, but the employed characteristic 
               particle weight is the average over all particles. 
          * 2: Bandwidth selection is performed taking into account only the particle 
               coordinates and a final density estimation is performed over the mass 
               histogram with the obtain distribution of kernel sizes.
          * 3: An effective number of particles is computed for each histogram cell and 
               is employed to transform the mass distribution into an effective particle 
               distribution used for the bandwidth optimization. A final reconstruction
               stage is performed considering the obtained distribution of kernel sizes.
    extension : str, optional
        File extension (default is 'gpkde').
    """

    
    @staticmethod
    def _ftype():
        return 'GPKDE'


    def __init__(
        self,
        model,
        outputfilename         = None,
        outputcolformat        = 0, 
        outputfileformat       = 0, 
        domainorigin           = None, 
        domainsize             = [1,1,1],
        binsize                = [1,1,1],
        gridallocformat        = 0, 
        gridborderfraction     = 0.05,
        noptloops              = 5,
        skiperror              = False, 
        convergence            = 0.02,
        kerneldatabase         = False, 
        isotropickernels       = False,
        boundkernelsize        = 0, 
        minhd                  = 1.0 ,
        maxhd                  = 0.1 ,
        deltahd                = 10.0,
        initialsmoothingformat = 0, 
        binsizefactor          = 5.0, 
        asconcentration        = False,
        effectiveweightformat  = 0,
        extension              = 'gpkde',
    ):


        ftype = self._ftype()
        unitnumber = model.next_unit()
        super().__init__(model, extension, ftype, unitnumber)
        
        # outputfilename
        if outputfilename is not None:
            if ( 
                not isinstance( outputfilename, str )  
            ):
                raise TypeError(
                    f"{self.__class__.__name__}:" 
                    f" The type of the output file name should be str. "
                    f"{str(type(outputfilename))} was given."
                )
            self.outputfilename = outputfilename
        else:
            self.outputfilename = f"{model.name}.gpkde.out"

        # outputcolformat
        if outputcolformat not in [0,1,2]: 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid value for outputcolformat. Allowed values are "
                f"0 (bin ids and density), 1 (bin ids, coordinates and density) "
                f"and 2 (coordinates and density). "
                f"{str(outputcolformat)} was given."
            )
        self.outputcolformat = outputcolformat

        # outfileformat
        if outputfileformat not in [0,1]: 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid value for outputfileformat. Allowed values are"
                f" 0 (text-plain) or 1 (binary)."
                f" {str(outputfileformat)}  was given."
            )
        self.outputfileformat = outputfileformat

        # domainorigin
        if domainorigin is None :
            # Fills domain origin with some values 
            # determined from the modelgrid definition
            domainorigin = [
                    self._parent.flowmodel.modelgrid.xoffset,
                    self._parent.flowmodel.modelgrid.yoffset,
                    np.min(self._parent.flowmodel.modelgrid.botm)
                ]
        self.domainorigin = np.array(domainorigin).astype(np.float32)
        if ( len(self.domainorigin.shape) != 1 ):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid specification of domain origin."
                f" It should be a one dimensional list with 3 elements."
            )
        if ( self.domainorigin.shape[0] != 3 ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid specification of domain origin."
                f" It should be a one dimensional list with 3 elements."
            )

        # Note: Depending on the kind of modelgrid maybe some
        # decisions can be made for binsizes and/or domain sizes.
        # StructuredGrid has the is_regular function.

        # domainsize
        if ( not isinstance( domainsize, (list,np.ndarray) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:" 
                f" The type of the domainsize parameter"
                f" should be list or np.ndarray."
                f" {str(type(domainsize))} was given."
            )
        self.domainsize = np.array(domainsize).astype(np.float32)
        if ( len(self.domainsize.shape) != 1 ):
            raise ValueError(
                f"{self.__class__.__name__}:" 
                f" Invalid specification of domain size."
                f" It should be a one dimensional list with 3 elements."
            )
        if ( self.domainsize.shape[0] != 3 ):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid specification of domain size."
                f" It should be a one dimensional list with 3 elements."
            )
        if ( np.any( self.domainsize < 0 ) ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Bin sizes should be greater than zero."
                f" {str(self.domainsize)} was given."
            )
        if ( np.all( self.domainsize == 0 ) ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" All domain sizes are zero."
                f" At least one positive value should be given."
            )

        # binsize
        if ( not isinstance( binsize, (list,np.ndarray) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" The type of the binsize parameter should be list or np.ndarray." 
                f" {str(type(binsize))}  was given."
            )
        self.binsize = np.array(binsize).astype(np.float32)
        if ( len(self.binsize.shape) != 1 ):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid specification of bin size."
                f" It should be a one dimensional list with 3 elements."
            )
        if ( self.binsize.shape[0] != 3 ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid specification of bin size."
                f" It should be a one dimensional list with 3 elements."
            )
        if ( np.any( self.binsize < 0 ) ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Bin sizes should be greater than zero."
                f" {str(self.binsize)} was given."
            )
        if ( np.all( self.binsize == 0 ) ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" All bin sizes are zero."
                f" At least one positive value should be given."
            )

        # gridallocformat
        if (gridallocformat not in [0,1]):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid value for gridallocformat. Allowed values are"
                f" 0 (allocate the grid as the domain) or 1 (allocate adapting to the particles)."
                f" {str(gridallocformat)} was given."
            )
        self.gridallocformat = gridallocformat

        # gridborderfraction
        if ( self.gridallocformat ==  1 ):
            # validate only if gridallocformat = 1
            if ( (not 0<gridborderfraction<=1) ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid value for gridborderfraction. Should be between 0 and 1."
                    f" 0 (allocate the grid as the domain) or 1 (allocate adapting to the particles)."
                    f" {str(gridborderfraction)} was given."
                )
            self.gridborderfraction = gridborderfraction

        # noptloops
        if ( not isinstance( noptloops, int) ): 
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for the number of optimization loops. It should be int, but"
                f" {str(type(noptloops))} was given."
            )
        if ( noptloops < 0 ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid value for the number of optimization loops. It should be a positive integer, but"
                f" {str(noptloops)} was given."
            )
        self.noptloops = noptloops

        # skiperror
        if ( not isinstance( skiperror, bool ) ): 
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for the skiperror parameter. It should be a boolean, but"
                f" {str(type(skiperror))} was given."
            )
        self.skiperror = skiperror

        # convergence
        if ( not isinstance( convergence, (int,float) ) ): 
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for the number of optimization loops. It should be float, but"
                f" {str(type(convergence))} was given."
            )
        if ( convergence < 0 ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid value for the convergence parameter. It should positive, but"
                f" {str(convergence)} was given."
            )
        self.convergence = convergence

        # kerneldatabase
        if ( not isinstance( kerneldatabase, bool ) ): 
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for the kerneldatabase parameter. It should be a boolean, but"
                f" {str(type(kerneldatabase))} was given."
            )
        self.kerneldatabase = kerneldatabase

        # isotropickernels
        if ( not isinstance( isotropickernels, bool ) ): 
            raise TypeError(
                f"{self.__class__.__name__}:" 
                f" Invalid type for the isotropickernels parameter. It should be a boolean, but"
                f" {str(type(isotropickernels))} was given."
            )
        self.isotropickernels = isotropickernels

        # boundkernelsize
        if ( boundkernelsize not in [0,1,2] ): 
            raise ValueError(
                f"{self.__class__.__name__}:" 
                f" Invalid value for boundkernelsize. Allowed values are 0 (bound by domain constraints),"
                f" 1 (bound by minhd,maxhd) or 2 (unbounded)."
                f" {str(boundkernelsize)} was given."
            )
        self.boundkernelsize = boundkernelsize

        # minhd
        # Verify only if required 
        if ( ( self.kerneldatabase ) or ( self.boundkernelsize==1 ) ):
            if ( not isinstance( minhd, (float,int) ) ):
                raise TypeError(
                    f"{self.__class__.__name__}:" 
                    f" Invalid type for minhd. It should be float or int, but"
                    f" {str(type(minhd))} was given."
                )
            if ( minhd <= 0 ):
                raise ValueError(
                    f"{self.__class__.__name__}:" 
                    f" Invalid value for minhd. It should be a positive value, but"
                    f" {str(type(minhd))} was given."
                )
            self.minhd = minhd

        # deltahd
        # Verify only if required 
        if ( ( self.kerneldatabase ) ):
            if ( not isinstance( deltahd, (float,int) ) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for deltahd. It should be float or int, but"
                    f" {str(type(deltahd))} was given."
                )
            if ( deltahd <= 0 ):
                raise ValueError(
                    f"{self.__class__.__name__}:" 
                    f" Invalid value for deltahd. It should be a positive value, but"
                    f" {str(type(deltahd))} was given."
                )
            self.deltahd = deltahd

        # maxhd
        # Verify only if required 
        if ( ( self.kerneldatabase ) or ( self.boundkernelsize==1 ) ):
            if ( not isinstance( maxhd, (float,int) ) ):
                raise TypeError(
                    f"{self.__class__.__name__}:" 
                    f" Invalid type for maxhd. It should be float or int, but"
                    f" {str(type(maxhd))} was given."
                )
            if ( maxhd <= 0 ):
                raise ValueError(
                    f"{self.__class__.__name__}:" 
                    f" Invalid value for maxhd. It should be a positive value, but"
                    f" {str(type(maxhd))} was given."
                )
            self.maxhd = maxhd

        # Minor consistency checks
        if ( ( self.kerneldatabase ) or ( self.boundkernelsize==1 ) ):
            if ( not( self.minhd <= self.maxhd ) ): 
                raise ValueError(
                    f"{self.__class__.__name__}:" 
                    f" Inconsistent values for minhd and maxhd."
                    f" Maxhd should be greater than minhd."
                )
            if ( self.kerneldatabase ):
                if ( ( self.deltahd > self.maxhd ) ): 
                    raise ValueError(
                        f"{self.__class__.__name__}:" 
                        f" Inconsistent values for deltahd and maxhd."
                        f" The given step value is higher than the maxhd."
                    )

        # initialsmoothingformat
        if (initialsmoothingformat not in [0,1]):
            raise ValueError(
                f"{self.__class__.__name__}:" 
                f" Invalid value for initialsmoothingformat. Allowed values are " 
                f" 0 (automatic initial selection) or 1 (initial kernel as a factor amplifying the bin size)."
                f" {str(initialsmoothingformat)} was given."
            )
        self.initialsmoothingformat = initialsmoothingformat

        # binsizefactor
        # Validate only if necessary
        if ( self.initialsmoothingformat == 1 ): 
            if ( not isinstance( binsizefactor, ( int, float ) ) ):
                raise TypeError(
                    f"{self.__class__.__name__}:" 
                    f" Invalid type for binsizefactor. It should be float or int, but"
                    f" {str(type(binsizefactor))} was given."
                )
            if ( binsizefactor <= 0 ):
                raise ValueError(
                    f"{self.__class__.__name__}:" 
                    f" Invalid value for binsizefactor. It should be a positive value, but"
                    f" {str(binsizefactor)} was given."
                )

        # asconcentration
        if ( not isinstance( asconcentration, bool ) ): 
            raise TypeError(
                f"{self.__class__.__name__}:" 
                f" Invalid type for the asconcentration parameter. It should be a boolean, but"
                f" {str(type(asconcentration))} was given."
            )
        self.asconcentration = asconcentration

        # effectiveweightformat
        if effectiveweightformat not in [0,1,2,3]: 
            raise ValueError(
                f"{self.__class__.__name__}:" 
                f" Invalid value for effectiveweightformat. This option determines the transformation to"
                f" an equivalent count of particles. Allowed values are 0 (unique effective mass),"
                f" 1 (scale by averaged mass), 2 (real histogram) or 3 (local effective histogram)."
                f" {str(effectiveweightformat)} was given."
            )
        self.effectiveweightformat = effectiveweightformat


        # Add package 
        self.parent.add_package(self)


        # And done
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

        # open file for writing
        with open(self.fn_path, "w") as f:

            # line 0:
            # Output file
            if ( ( self.outputcolformat != 0 ) or ( self.outputfileformat != 0 ) ): 
                f.write(f"{self.outputfilename}  {self.outputcolformat}  {self.outputfileformat}\n")
            else:
                f.write(f"{self.outputfilename}\n")

            # line 1:
            # domain origin
            for idb, b in enumerate(self.domainorigin):
                if idb == len(self.domainorigin)-1: 
                    f.write(f"{b:10f}\n")
                else:
                    f.write(f"{b:10f}\t")

            # line 2:
            # domain size
            if ( self.gridallocformat == 1 ):
                if ( self.gridborderfraction > 0 ):
                    for idb, b in enumerate(self.domainsize):
                        f.write(f"{b:10f}\t")
                    f.write(f"{self.gridallocformat}\t")
                    f.write(f"{self.gridborderfraction}\n")
                else:
                    for idb, b in enumerate(self.domainsize):
                        f.write(f"{b:10f}\t")
                    f.write(f"{self.gridallocformat}\n")
            else:
                for idb, b in enumerate(self.domainsize):
                    if idb == len(self.domainsize)-1: 
                        f.write(f"{b:10f}\n")
                    else:
                        f.write(f"{b:10f}\t")

            # line 3:
            # bin sizes
            for idb, b in enumerate(self.binsize):
                if idb == len(self.binsize)-1: 
                    f.write(f"{b:10f}\n")
                else:
                    f.write(f"{b:10f}\t")

            # line 4:
            # number of optimization loops 
            f.write(f"{self.noptloops}\n")

            # line 5:
            # skip error convergence
            if ( self.skiperror ): 
                f.write(f"{int(self.skiperror)}\n")
            else:
                f.write(f"{int(self.skiperror)}    {self.convergence:6f}\n")

            # line 6:
            # kernels
            if ( not self.isotropickernels ): 
                f.write(f"{int(self.kerneldatabase)}    {self.boundkernelsize}\n")
            else:
                f.write(f"{int(self.kerneldatabase)}    {self.boundkernelsize}    {int(self.isotropickernels)}\n")

            # line 7:
            # kdb params
            if ( self.kerneldatabase ):
                # Database params: minhd, deltahd, maxhd
                f.write(f"{self.minhd:10f}\t")
                f.write(f"{self.deltahd:10f}\t")
                f.write(f"{self.maxhd:10f}\n")
            elif ( self.boundkernelsize == 1):
                # minhd, maxhd
                f.write(f"{self.minhd:10f}\t")
                f.write(f"{self.maxhd:10f}\n")

            # line 8:
            # initial smoothing
            if ( self.initialsmoothingformat == 1 ):
                f.write(f"{self.initialsmoothingformat}   {self.binsizefactor}\n") 
            else:
                f.write(f"{self.initialsmoothingformat}\n") 

            # line 9
            # as concentration 
            f.write(f"{int(self.asconcentration)}\n")

            # line 10
            # as concentration 
            f.write(f"{self.effectiveweightformat}\n")

        # Done
        return
