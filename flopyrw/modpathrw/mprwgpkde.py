'''
Configuration of MODPATH-RW gpkde reconstruction
'''

# python
import numpy as np
from enum import Enum

# flopy
from flopy.utils import Util2d
from flopy.pakbase import Package
from flopy.discretization import StructuredGrid, VertexGrid, UnstructuredGrid


class slicedDimensionOption(Enum): 
    '''
    Enumerate possible values for sliced dimension 
    '''
    x = 0
    y = 1
    z = 2


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
    slicedreconstruction : bool
        Performs reconstruction for each slice of a given dimension. Recommended 
        for quasi two-dimensional problems, for example, when the horizontal extent
        dominates over the vertical. 
    sliceddimension: int or str
        The dimension to be sliced. Allowed values are:
          0 or 'x': reconstruction for each slice of dimension x 
          1 or 'y': reconstruction for each slice of dimension y 
          2 or 'z': reconstruction for each slice of dimension z
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
    anisotropickernelsic : int
        Indicate the isotropy status for kernels employed in the reconstruction of the 
        initial condition. Allowed values are:
          * 0: kernels are isotropic.
          * 1: kernels are anisotropic.
          * 2: kernels follow the variable isotropickernels.
    boundkernelsize : int
        Determines the protocol for bounding the kernel size. Allowed values are: 
          * 0: bound size by domain constraints. 
          * 1: bound using minhd and maxhd
          * 2: unbounded
    minhd : float
        Minimum ratio smoothing(h)/binsize(l). It is written to the package file 
        if kerneldatabase=True or boundkernelsize = 1.
    deltahd : float
        Step of ratio smoothing(h)/binsize(l) used for defining a set of discrete
        kernel sizes for the kernel database. Only written to the package file 
        if kerneldatabase=True.
    maxhd : float
        Maximum ratio smoothing(h)/binsize(l). It is written to the package file if
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
    timepointdata : list or tuple
        Follows a similar logic to timepointdata at the ModpathRWSim/Modpath7Sim class.
        Is a list or tuple with 2 items. If the second item is a float then 
        thetimepointdata corresponds to timepointoption 1 and the first entry is the
        number of time points (timepointcount) and the second the interval.
        If the second item is a list, tuple, or np.ndarray, timepointoption is 2 and
        reconstruction is performed for the times of the given list. If timepointdata is None,
        timepointoption is 0 and reconstruction follows the timeseries points. 
    skipinitialcondition: bool
        Flag to indicate whether reconstruction should be performed for the 
        initial distribution of particles.
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
        domainsize             = None, 
        binsize                = None,
        gridallocformat        = 0, 
        gridborderfraction     = 0.05,
        slicedreconstruction   = None, 
        sliceddimension        = None, 
        noptloops              = 10,
        skiperror              = False, 
        convergence            = 0.02,
        kerneldatabase         = False, 
        isotropickernels       = False,
        boundkernelsize        = 0, 
        anisotropickernelsic   = 0, 
        minhd                  = 1.0 ,
        maxhd                  = 0.1 ,
        deltahd                = 10.0,
        initialsmoothingformat = 0, 
        binsizefactor          = 1.0, 
        asconcentration        = True,
        effectiveweightformat  = 0,
        timepointdata          = None,
        skipinitialcondition   = False,
        extension              = 'gpkde',
    ):
        
        # initialize
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
                    f" The type of the output file name should be str."
                    f" {str(type(outputfilename))} was given."
                )
            self.outputfilename = outputfilename
        else:
            self.outputfilename = f"{model.name}.gpkde.out"

        # outputcolformat
        if outputcolformat not in [0,1,2]: 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid value for outputcolformat. Allowed values are"
                f" 0 (bin ids and density), 1 (bin ids, coordinates and density)"
                f" and 2 (coordinates and density)."
                f" {str(outputcolformat)} was given."
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


        # Get the modelgrid for convenience
        modelgrid = self.parent.flowmodel.modelgrid 

        # domainsize
        if (
            ( domainsize is None ) and 
            isinstance( modelgrid, (StructuredGrid,VertexGrid) )
        ):

            domainsize = np.zeros(shape=(3,)).astype(np.float32)
            xmin,xmax,ymin,ymax = modelgrid.extent

            if ( isinstance( modelgrid, StructuredGrid ) ):
                # Assign domainsizes
                domainsize[0] = abs( xmax - xmin )
                domainsize[1] = abs( ymax - ymin )
                zmin = np.min( modelgrid.botm )
                zmax = np.max( modelgrid.top  ) 
                domainsize[2] = abs( zmax - zmin )
            elif ( isinstance( modelgrid, VertexGrid ) ):
                # Is the same than for structured 
                domainsize[0] = abs( xmax - xmin )
                domainsize[1] = abs( ymax - ymin )
                zmin = np.min( modelgrid.botm )
                zmax = np.max( modelgrid.top  ) 
                domainsize[2] = abs( zmax - zmin )
        elif (
            ( domainsize is None ) and 
            isinstance( modelgrid, (UnstructuredGrid) )
        ):
            # In the meantime, force explicit specification
            if ( self.parent.flowmodel.modelversion == 'mf6' ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" For unstructured grids an explicit specification"
                    f" of domainsize is needed, but None was given."
                    f" Note: Modpath7 and ModpathRW do not currently support"
                    f" unstructured grids for mf6 flow models."
                )
            else: 
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" For unstructured grids an explicit specification"
                    f" of domainsize is needed, but None was given."
                )

        if ( not isinstance( domainsize, (list,np.ndarray) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:" 
                f" The type of the domainsize parameter"
                f" should be list or np.ndarray."
                f" {str(type(domainsize))} was given."
            )
        domainsize = np.array(domainsize).astype(np.float32)
        if ( len(domainsize.shape) != 1 ):
            raise ValueError(
                f"{self.__class__.__name__}:" 
                f" Invalid specification of domain size."
                f" It should be a one dimensional list with 3 elements."
            )
        if ( domainsize.shape[0] != 3 ):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid specification of domain size."
                f" It should be a one dimensional list with 3 elements."
            )
        if ( np.any( domainsize < 0 ) ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Bin sizes should be greater than zero."
                f" {str(self.domainsize)} was given."
            )
        if ( np.all( domainsize == 0 ) ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" All domain sizes are zero."
                f" At least one positive value should be given."
            )
        self.domainsize = domainsize


        # domainorigin
        if (
            ( domainorigin is None ) and 
            isinstance( modelgrid, (StructuredGrid,VertexGrid) )
        ):

            domainorigin = np.zeros(shape=(3,)).astype(np.float32)

            if ( isinstance( modelgrid, StructuredGrid ) ):
                # Fills domain origin with values 
                # determined from the modelgrid definition
                domainorigin[0] = modelgrid.xoffset
                domainorigin[1] = modelgrid.yoffset
                domainorigin[2] = np.min(modelgrid.botm)
            elif ( isinstance( modelgrid, VertexGrid ) ):
                # The same as structured 
                domainorigin[0] = modelgrid.xoffset
                domainorigin[1] = modelgrid.yoffset
                domainorigin[2] = np.min(modelgrid.botm)
        elif (
            ( domainorigin is None ) and 
            isinstance( modelgrid, (UnstructuredGrid) )
        ):
            # In the meantime, force explicit specification
            if ( self.parent.flowmodel.modelversion == 'mf6' ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" For unstructured grids an explicit specification"
                    f" of domainorigin is needed, but None was given."
                    f" Note: Modpath7 and ModpathRW do not currently support"
                    f" unstructured grids for mf6 flow models."
                )
            else: 
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" For unstructured grids an explicit specification"
                    f" of domainorigin is needed, but None was given."
                )

        if ( not isinstance( domainorigin, (list,np.ndarray) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:" 
                f" The type of the domainorigin parameter"
                f" should be list or np.ndarray."
                f" {str(type(domainorigin))} was given."
            )
        domainorigin = np.array(domainorigin).astype(np.float32)
        if ( len(domainorigin.shape) != 1 ):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid specification of domain origin."
                f" It should be a one dimensional list with 3 elements."
            )
        if ( domainorigin.shape[0] != 3 ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid specification of domain origin."
                f" It should be a one dimensional list with 3 elements."
            )
        self.domainorigin = domainorigin


        # binsize
        if (
            ( binsize is None ) and 
            isinstance( modelgrid, (StructuredGrid,VertexGrid) )
        ): 
            # For a structured grid, infer bin sizes
            # from delr,delc,delz
            binsize = np.zeros(shape=(3,)).astype(np.float32)

            if ( isinstance( modelgrid, StructuredGrid ) ):
                # x
                if modelgrid.is_regular_x:
                    # If the grid is regular, then get the value
                    binsize[0] = modelgrid.delr.flat[0]
                else:
                    # If not, it could be an average
                    binsize[0] = np.mean(modelgrid.delr.flat)
                # y
                if modelgrid.is_regular_y:
                    binsize[1] = modelgrid.delc.flat[0]
                else:
                    binsize[1] = np.mean(modelgrid.delc.flat)
                # z
                if modelgrid.is_regular_z:
                    binsize[2] = modelgrid.delz.flat[0]
                else:
                    binsize[2] = np.mean(modelgrid.delz.flat)
            elif ( isinstance( modelgrid, VertexGrid ) ):
                # Infer a bin size based on xcellcenters 
                # and ycellcenters
                # x
                # Check if unique difference
                xunique = np.unique(modelgrid.xcellcenters)
                if len(xunique) > 1:
                    xdiff = np.diff(xunique)
                    if ( len(np.unique(xdiff))==1 ):
                        # Take as bin size if unique
                        binsize[0] = np.unique(xdiff).item()
                    else:
                        # otherwise an average
                        binsize[0] = np.mean(np.unique(xdiff))
                else:
                    # if one element, assign domainsize
                    binsize[0] = domainsize[0]
                # y
                yunique = np.unique(modelgrid.ycellcenters)
                if len(yunique) > 1:
                    ydiff = np.diff(yunique)
                    if ( len(np.unique(ydiff))==1 ):
                        binsize[1] = np.unique(ydiff).item()
                    else:
                        binsize[1] = np.mean(np.unique(ydiff))
                else:
                    # if one, assign domainsize
                    binsize[1] = domainsize[1]
                # z
                zunique = np.unique(modelgrid.zcellcenters)
                if len(zunique) > 1:
                    zdiff = np.diff(zunique)
                    if ( len(np.unique(zdiff))==1 ):
                        binsize[2] = np.unique(zdiff).item()
                    else:
                        binsize[2] = np.mean(np.unique(zdiff))
                else:
                    # if one, assign domainsize
                    binsize[2] = domainsize[2]
        elif (
            ( binsize is None ) and 
            isinstance( modelgrid, (UnstructuredGrid) )
        ):
            # In the meantime, force explicit specification
            if ( self.parent.flowmodel.modelversion == 'mf6' ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" For unstructured grids an explicit specification"
                    f" of binsize is needed, but None was given."
                    f" Note: Modpath7 and ModpathRW do not currently support"
                    f" unstructured grids for mf6 flow models."
                )
            else: 
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" For unstructured grids an explicit specification"
                    f" of binsize is needed, but None was given."
                )

        if ( not isinstance( binsize, (list,np.ndarray) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" The type of the binsize parameter should be list or np.ndarray." 
                f" {str(type(binsize))}  was given."
            )
        binsize = np.array(binsize).astype(np.float32)
        if ( len(binsize.shape) != 1 ):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid specification of bin size."
                f" It should be a one dimensional list with 3 elements."
            )
        if ( binsize.shape[0] != 3 ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid specification of bin size."
                f" It should be a one dimensional list with 3 elements."
            )
        if ( np.any( binsize < 0 ) ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Bin sizes should be greater than zero."
                f" {str(self.binsize)} was given."
            )
        if ( np.all( binsize == 0 ) ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" All bin sizes are zero."
                f" At least one positive value should be given."
            )
        self.binsize = binsize

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

        # sliced reconstruction
        if ( slicedreconstruction is None ): 
            # some inference for 3d problems
            if ( np.all( binsize >  0 ) ):
                # rough estimate of the number of bins
                nbins = domainsize/binsize
                # skip inference if any nbins is one, 
                # gpkde will recognize this case.
                if np.any( nbins == 1 ): 
                    self.slicedreconstruction = False
                    self.sliceddimension = 0
                else:
                    # Get the minimum number of bins and estimate 
                    # anisotropy
                    idmin, = np.where( nbins == np.min( nbins ) )
                    idmin  = idmin.item()
                    # none of these seems to be that general
                    ani13  = nbins[(idmin+1)%3]/nbins[idmin]
                    ani23  = nbins[(idmin+2)%3]/nbins[idmin]
                    #ani13  = domainsize[(idmin+1)%3]/domainsize[idmin]
                    #ani23  = domainsize[(idmin+2)%3]/domainsize[idmin]

                    # Define an arbitrary threshold of 10
                    if ( ( ani13 > 10 ) and ( ani23 > 10 ) ) :
                        self.slicedreconstruction = True
                        self.sliceddimension = idmin
                    else:
                        # the inference might be improved, but in the meantime
                        # continue as false.
                        self.slicedreconstruction = False
                        self.sliceddimension = 0
            else:
                self.slicedreconstruction = False
                self.sliceddimension = 0
        else:
            if ( not isinstance( slicedreconstruction, bool ) ) : 
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for the slicedreconstruction parameter. It should be bool, but"
                    f" {str(type(slicedreconstruction))} was given."
                )
            self.slicedreconstruction = slicedreconstruction

            # process sliceddimension
            if ( self.slicedreconstruction ): 
                if ( sliceddimension is None ):
                    # Do the inference as before
                    # some inference for 3d problems
                    if ( np.all( binsize >  0 ) ):
                        # rough estimate of the number of bins
                        nbins = domainsize/binsize
                        # skip inference if any nbins is one, 
                        # gpkde will recognize this case.
                        if np.any( nbins == 1 ):
                            # force the user to provide a value
                            raise ValueError(
                                f"{self.__class__.__name__}:"
                                f" Invalid value for sliceddimension. None was provided and"
                                f" inference could not be done."
                            )

                        else:
                            # Get the minimum number of bins and estimate 
                            # anisotropy
                            idmin, = np.where( nbins == np.min( nbins ) )
                            idmin  = idmin.item()
                            # none of these seems to be that general
                            ani13  = nbins[(idmin+1)%3]/nbins[idmin]
                            ani23  = nbins[(idmin+2)%3]/nbins[idmin]
                            #ani13  = domainsize[(idmin+1)%3]/domainsize[idmin]
                            #ani23  = domainsize[(idmin+2)%3]/domainsize[idmin]

                            # Define an arbitrary threshold of 10
                            if ( ( ani13 > 10 ) and ( ani23 > 10 ) ) :
                                self.slicedreconstruction = True
                                self.sliceddimension = idmin
                            else:
                                # the inference might be improved, but in the meantime
                                # force the user to provide a value
                                raise ValueError(
                                    f"{self.__class__.__name__}:"
                                    f" Invalid value for sliceddimension. None was provided and"
                                    f" inference could not be done."
                                )

                elif ( not isinstance( sliceddimension, (int,str) ) ):
                    raise TypeError(
                        f"{self.__class__.__name__}:"
                        f" Invalid type for the sliceddimension parameter. It should be int or str, but"
                        f" {str(type(sliceddimension))} was given."
                    )
                else:
                    if isinstance( sliceddimension, int ) :
                        if ( sliceddimension not in [0,1,2] ):
                            raise ValueError(
                                f"{self.__class__.__name__}:"
                                f" Invalid sliceddimension {str(sliceddimension)}."
                                f" The allowed values are 0 (x), 1 (y) or 2 (z)."
                            )
                        self.sliceddimension = sliceddimension + 1 
                    else:
                        try:
                            self.sliceddimension = slicedDimensionOption[sliceddimension.lower()].value + 1
                        except:
                            raise ValueError(
                                f"{self.__class__.__name__}:"
                                f" Invalid sliceddimension {str(sliceddimension)}."
                                f" The allowed values are x, y or z."
                            )

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

        # anisotropickernelsic
        if ( not isinstance( anisotropickernelsic, int ) ): 
            raise TypeError(
                f"{self.__class__.__name__}:" 
                f" Invalid type for the anisotropickernelsic parameter. It should be int, but"
                f" {str(type(anisotropickernelsic))} was given."
            )
        if ( anisotropickernelsic not in [0,1,2] ): 
            raise ValueError(
                f"{self.__class__.__name__}:" 
                f" Invalid value for anisotropickernelsic . Allowed values are 0 (isotropic kernels for ic),"
                f" 1 (anisotropic kernels for ic) or 2 (follow isotropickernels)."
                f" {str(anisotropickernelsic)} was given."
            )
        self.anisotropickernelsic = anisotropickernelsic

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
            self.binsizefactor = binsizefactor

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

        # timepointdata, similar logic to mp7sim
        if timepointdata is not None:
            if not isinstance(timepointdata, (list, tuple)):
                raise TypeError(
                        f"{self.__class__.__name__}:" 
                        f" Invalid type for timepointdata, it"
                        f" should be a list or tuple, but"
                        f" {type(timepoindata)} was given."
                    )
            else:
                if len(timepointdata) != 2:
                    raise ValueError(
                        f"{self.__class__.__name__}:" 
                        f" timepointdata must be a have 2 entries, but "
                        f" {len(timepointdata)} were provided."
                    )
                else:
                    # Infers the specification
                    if isinstance(timepointdata[1], (list, tuple)):
                        timepointdata[1] = np.array(timepointdata[1])
                    elif isinstance(timepointdata[1], (int,float)):
                        timepointdata[1] = np.array([timepointdata[1]]).astype(np.float32)

                    # Infers timepointoption
                    if timepointdata[1].shape[0] == timepointdata[0]:
                        timepointoption = 2
                    elif timepointdata[1].shape[0] > 1:
                        raise ValueError(
                            f"{self.__class__.__name__}:" 
                            f" The number of given timepoints {str(timepointdata[1].shape[0])}"
                            f" is not equal to the given timepointcount {str(timepointdata[0])}"
                        )
                    else:
                        timepointoption = 1
        else:
            # follows timeseries
            timepointoption = 0
        self.timepointoption = timepointoption
        self.timepointdata   = timepointdata

        # skipinitialcondition
        if ( not isinstance( skipinitialcondition, bool ) ): 
            raise TypeError(
                    f"{self.__class__.__name__}:" 
                    f" Invalid type for skipinitialcondition, it"
                    f" should be boolean, but "
                    f" {type(skipinitialcondition)} was given."
                )
        if skipinitialcondition: 
            self.skipinitialcondition = 1
        else:
            self.skipinitialcondition = 0

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
            # sliced reconstruction 
            if not self.slicedreconstruction:
                f.write(f"0\n")
            else:
                f.write(f"1      {self.sliceddimension}\n")

            # line 5:
            # number of optimization loops 
            f.write(f"{self.noptloops}\n")

            # line 6:
            # skip error convergence
            if ( self.skiperror ): 
                f.write(f"{int(self.skiperror)}\n")
            else:
                f.write(f"{int(self.skiperror)}   {self.convergence:6f}\n")

            # line 7:
            # kernels
            f.write(f"{int(self.kerneldatabase)}   {self.boundkernelsize}   {int(self.isotropickernels)}   {int(self.anisotropickernelsic)}\n")

            # line 8:
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

            # line 9:
            # initial smoothing
            if ( self.initialsmoothingformat == 1 ):
                f.write(f"{self.initialsmoothingformat}   {self.binsizefactor}\n") 
            else:
                f.write(f"{self.initialsmoothingformat}\n") 

            # line 10
            # as concentration 
            f.write(f"{int(self.asconcentration)}\n")

            # line 11
            # effectiveweightformat
            f.write(f"{self.effectiveweightformat}\n")

            # line 12 
            # timepointoption and skipinitialcondition
            f.write(f"{self.timepointoption}    {self.skipinitialcondition}\n")
           
            # timepointdata
            if self.timepointoption == 1:
                # line 13
                f.write(f"{self.timepointdata[0]}    {self.timepointdata[1][0]}\n")
            elif self.timepointoption == 2:
                # line 13, 14
                f.write(f"{self.timepointdata[0]}\n")
                tp = self.timepointdata[1]
                v  = Util2d(
                    self.parent,
                    (tp.shape[0],),
                    np.float32,
                    tp,
                    name="tpoints",
                    locat=self.unit_number[0],
                )
                f.write(v.string)

        # Done
        return


    def get_output(self):
        '''
        Load the output file into a rec array
        '''
        import os 
        from flopy.utils.flopy_io import loadtxt
        # Shouldn't this be a case to handle by loadtxt ?
        try:
            import pandas as pd
            use_pandas = True
        except ImportError:
            use_pandas = False
        #
        # -- define dtype and vartype 
        #    depending on the output column format
        if self.outputcolformat == 0:
            dtype = np.dtype(
                [
                    ("tid"       , np.int32  ),
                    ("time"      , np.float32),
                    ("speciesid" , np.int32  ),
                    ("idbinx"    , np.int32  ),
                    ("idbiny"    , np.int32  ),
                    ("idbinz"    , np.int32  ),
                    ("cgpkde"    , np.float32),
                    ("chist"     , np.float32),
                ]
            )
            vartype = [
                ("tid"      , "<i4"),
                ("time"     , "<f8"),
                ("speciesid", "<i4"),
                ("idbinx"   , "<i4"),
                ("idbiny"   , "<i4"),
                ("idbinz"   , "<i4"),
                ("cgpkde"   , "<f8"),
                ("chist"    , "<f8"),
            ]
        elif self.outputcolformat == 1:
            dtype = np.dtype(
                [
                    ("tid"       , np.int32  ),
                    ("time"      , np.float32),
                    ("speciesid" , np.int32  ),
                    ("idbinx"    , np.int32),
                    ("idbiny"    , np.int32),
                    ("idbinz"    , np.int32),
                    ("x"         , np.float32),
                    ("y"         , np.float32),
                    ("z"         , np.float32),
                    ("cgpkde"    , np.float32),
                    ("chist"     , np.float32),
                ]
            )
            vartype = [
                ("tid"      , "<i4"),
                ("time"     , "<f8"),
                ("speciesid", "<i4"),
                ("idbinx"   , "<i4"),
                ("idbiny"   , "<i4"),
                ("idbinz"   , "<i4"),
                ("x"        , "<f8"),
                ("y"        , "<f8"),
                ("z"        , "<f8"),
                ("cgpkde"   , "<f8"),
                ("chist"    , "<f8"),
            ]
        elif self.outputcolformat == 2:
            dtype = np.dtype(
                [
                    ("tid"       , np.int32  ),
                    ("time"      , np.float32),
                    ("speciesid" , np.int32  ),
                    ("x"         , np.float32),
                    ("y"         , np.float32),
                    ("z"         , np.float32),
                    ("cgpkde"    , np.float32),
                    ("chist"     , np.float32),
                ]
            )
            vartype = [
                ("tid"      , "<i4"),
                ("time"     , "<f8"),
                ("speciesid", "<i4"),
                ("x"        , "<f8"),
                ("y"        , "<f8"),
                ("z"        , "<f8"),
                ("cgpkde"   , "<f8"),
                ("chist"    , "<f8"),
            ]
        #
        # -- read data
        if self.outputfileformat == 0:
            # -- from text-plain file
            recdata = loadtxt(
                os.path.join( self._parent.model_ws, self.outputfilename ),
                dtype=dtype,
                skiprows=0, 
                use_pandas=use_pandas,
            )
        elif self.outputfileformat == 1:
            # -- from binary file
            file = open( os.path.join( self._parent.model_ws, self.outputfilename ), 'rb' )
            data = np.fromfile(file, vartype)
            recdata = data.view(np.recarray)
        #
        # -- to python zero-based indexes those quantities 
        #    that need this treatment. leave tid as it came.
        recdata['speciesid'] = recdata['speciesid'] - 1
        if self.outputcolformat != 2:
            recdata['idbinx'] = recdata['idbinx'] - 1
            recdata['idbiny'] = recdata['idbiny'] - 1
            recdata['idbinz'] = recdata['idbinz'] - 1
        # 
        # -- store variables
        self.times = np.unique(recdata['time'])
        self.speciesids = np.unique(recdata['speciesid'])
        #
        # -- useful, heavy ?
        self.outputrecdata = recdata
        #
        # -- return
        return recdata 


    def get_times(self):
        '''
        Return the array of times stored in data file
        '''
        try:
            return self.times
        except AttributeError:
            self.get_output()
            return self.times


    def get_data(self, which='cgpkde', totim=None, speciesid=None):
        '''
        Get output data array

        For a reconstruction grid coincident with the flowmodel grid of 
        type StructuredGrid and regular cell sizes, will fill an array with shape 
        (nlay,nrow,ncol), with the concentration data requested in 'which'.
          * which=cgpkde returns the smoothed reconstruction. 
          * which=chist returns the histogram reconstruction. 

        In case the grid is not of type StructuredGrid or not regular, 
        will filter the recarray data by speciesid and totim, and not by 'which'.
        The same applies for the case in which self.outputcolformat == 2, where 
        reconstruction grid indexes are not given in the output file. 

        Note: the gpkde reconstruction grid follows the modpath convention of 
              coordinates and not the lay,row,col convention so there should an
              adequate reinterpration of grid indexes.
        '''
    
        # Default value for non existent indexes
        defaultnan = 0.0

        # Validate which
        if not isinstance(which,str):
            raise TypeError(
                f"{self.__class__.__name__}:get_data:" 
                f" Invalid type for which parameter, it should be str but {str(type(which))} was given."
            )
        if (which.lower() not in ['cgpkde', 'chist']):
            raise ValueError(
                f"{self.__class__.__name__}:get_data:"
                f" Invalid value for which parameter. It can be"
                f" cgpkde or chist, but {str(which)} was given."
            )
        which = which.lower()

        # Load the output file if not loaded already
        try: 
            recdata = self.outputrecdata
        except AttributeError:
            recdata = self.get_output()

        # If total time was none, assign the last in 
        # self.times
        if totim is None: 
            totim = self.times[-1]
        else:
            tindex = np.where(self.times==totim)[0]
            if len(tindex)==0:
                raise ValueError(
                    f"{self.__class__.__name__}:get_data:"
                    f" The given value for totim was not found in the array of times. "
                    f" totim={str(totim)} was given."
                )
            else:
                totim = self.times[tindex.item()]

        # Similar to times, do it for speciesid
        if speciesid is None: 
            speciesid = self.speciesids[-1]
        else:
            spcindex = np.where(self.speciesids==speciesid)[0]
            if len(spcindex)==0:
                raise ValueError(
                    f"{self.__class__.__name__}:get_data:"
                    f" The given value for speciesid was not found in the array of speciesids. "
                    f" speciesid={str(speciesid)} was given."
                )
            else:
                speciesid = self.speciesids[spcindex.item()]

        # Filter data
        filtdata = recdata[ (recdata['time'] == totim)&(recdata['speciesid'] == speciesid) ]

        # If no data, error.
        if ( len(filtdata) == 0 ):
            raise Exception(
                f"{self.__class__.__name__}:get_data:"
                f" No data was found for totim={str(totim)} and speciesid={str(speciesid)}"
            )

        # If no grid indexes, return
        if ( self.outputcolformat == 2 ):
            return filtdata

        # 
        # The following would work for StructuredGrid with regular cell sizes. 
        # Additional alternatives could be provided, for example interpolating 
        # with griddata for unstructred grids.
        if isinstance( self.parent.flowmodel.modelgrid, StructuredGrid ):
            #
            # -- flag to determine whether it is possible to load 
            #    gpkde data with the same structure than a modflow 
            #    grid.
            canloaddata = False
            # 
            # -- evaluate the condition 
            if self.parent.flowmodel.modelgrid.is_regular:
                # -- if the grid is regular, cubic cells, then 
                #    it is possible 
                canloaddata = True
            elif (
                self.parent.flowmodel.modelgrid.is_regular_x and 
                self.parent.flowmodel.modelgrid.is_regular_y and 
                self.parent.flowmodel.modelgrid.is_regular_z  ):
                # -- if the grid is regular in each axis.
                #    this case is possible when cell sizes are 
                #    are not the same, but regular for each axis.
                canloaddata = True
            # 
            # -- here let's do a simple evaluation on the 
            #    number of cells
            nbins = np.zeros(shape=(3,))
            for ib in range(3):
                if self.binsize[ib] == 0:
                    continue
                # -- follows the x,y,z convention
                nbins[ib] = int(self.domainsize[ib]/self.binsize[ib])
            # 
            # -- get the modflow grid size
            nlay = self.parent.flowmodel.modelgrid.nlay
            nrow = self.parent.flowmodel.modelgrid.nrow
            ncol = self.parent.flowmodel.modelgrid.ncol
            # 
            # -- validate same grid size.
            #    if binsize is zero, then allow
            #    loading if the modflow grid has only
            #    one element in the given dimension
            # 
            # -- x 
            if (nbins[0] != 0):
                if ( nbins[0] != ncol ): 
                    canloaddata = False
            else:
                if ( ncol != 1 ): 
                    canloaddata = False
            # 
            # -- y
            if (nbins[1] != 0):
                if ( nbins[1] != nrow ): 
                    canloaddata = False
            else:
                if ( nrow != 1 ): 
                    canloaddata = False
            # 
            # -- z
            if (nbins[2] != 0):
                if ( nbins[2] != nlay ): 
                    canloaddata = False
            else:
                if ( nlay != 1 ): 
                    canloaddata = False
            #
            # -- load the gpkde data into an array
            if canloaddata:
                # 
                data = np.empty((nlay, nrow, ncol), dtype=np.float32)
                data[:,:,:] = defaultnan
                #
                # -- load for the requested time 
                layers = np.unique( filtdata['idbinz'] )
                for lay in layers:
                    srec = filtdata[ filtdata['idbinz'] == lay ]
                    laydata = np.empty((nrow,ncol),dtype=np.float32)
                    laydata[:,:] = defaultnan
                    laydata[ nrow - srec['idbiny'] - 1, srec['idbinx'] ] = srec[which] 
                    data[nlay - lay - 1,:,:] = laydata
                # -- return
                return data
            else:
                print( 
                    f"Warning: data for flowmodel with non-regular StructuredGrid" 
                    f" is returned filtered only by totim and speciesid. The user should apply"
                    f" an adequate coordinates conversion (e.g., scipy.interpolate.griddata)."
                )
                # -- return
                return filtdata
        else:
            print( 
                f"Warning: data for flowmodel with {str(type(self.parent.flowmodel.modelgrid))} grid"
                f" is returned filtered only by totim and speciesid. The user should apply"
                f" an adequate coordinates conversion. "
            )
            # return
            return filtdata


    def get_alldata(self, which='cgpkde', speciesid=None):
        '''
        Get output data array for all times

        For a reconstruction grid coincident with the flowmodel grid of 
        type StructuredGrid and regular cell sizes, will fill an array with shape 
        (ntimes,nlay,nrow,ncol), with the concentration data requested in 'which'.
          * which=cgpkde returns the smoothed reconstruction. 
          * which=chist returns the histogram reconstruction. 

        In case the grid is not of type StructuredGrid or not regular, 
        will filter the recarray data by speciesid, and not by 'which'.
        The same applies for the case in which self.outputcolformat == 2, where 
        reconstruction grid indexes are not given in the output file. 

        Note: the gpkde reconstruction grid follows the modpath convention of 
              coordinates and not the lay,row,col convention so there should an
              adequate reinterpration of grid indexes.
        '''
    
        # Default value for non existent indexes
        defaultnan = 0.0

        # Validate which
        if not isinstance(which,str):
            raise TypeError(
                f"{self.__class__.__name__}:get_alldata:" 
                f" Invalid type for which parameter, it should be str but {str(type(which))} was given."
            )
        if (which.lower() not in ['cgpkde', 'chist']):
            raise ValueError(
                f"{self.__class__.__name__}:get_alldata:"
                f" Invalid value for which parameter. It can be"
                f" cgpkde or chist, but {str(which)} was given."
            )
        which = which.lower()

        # Load the output file if not loaded already
        try: 
            recdata = self.outputrecdata
        except AttributeError:
            recdata = self.get_output()

        # get total number of times 
        ntimes = self.times.shape[0]

        # Get a speciesid, by default the last
        if speciesid is None: 
            speciesid = self.speciesids[-1]
        else:
            spcindex = np.where(self.speciesids==speciesid)[0]
            if len(spcindex)==0:
                raise ValueError(
                    f"{self.__class__.__name__}:get_alldata:"
                    f" The given value for speciesid was not found in the array of speciesids. "
                    f" speciesid={str(speciesid)} was given."
                )
            else:
                speciesid = self.speciesids[spcindex.item()]

        # Filter data: only by speciesid
        filtdata = recdata[ (recdata['speciesid'] == speciesid) ]

        # If no data, error.
        if ( len(filtdata) == 0 ):
            raise Exception(
                f"{self.__class__.__name__}:get_alldata:"
                f" No data was found for speciesid={str(speciesid)}"
            )
        
        # If no grid indexes, return
        if ( self.outputcolformat == 2 ):
            return filtdata
        # 
        # The following would work for StructuredGrid with regular cell sizes. 
        # Additional alternatives could be provided, for example interpolating 
        # with griddata for unstructred grids.
        if isinstance( self.parent.flowmodel.modelgrid, StructuredGrid ):
            #
            # -- flag to determine whether it is possible to load 
            #    gpkde data with the same structure than a modflow 
            #    grid.
            canloaddata = False
            # 
            # -- evaluate the condition 
            if self.parent.flowmodel.modelgrid.is_regular:
                # -- if the grid is regular, cubic cells, then 
                #    it is possible 
                canloaddata = True
            elif (
                self.parent.flowmodel.modelgrid.is_regular_x and 
                self.parent.flowmodel.modelgrid.is_regular_y and 
                self.parent.flowmodel.modelgrid.is_regular_z  ):
                # -- if the grid is regular in each axis.
                #    this case is possible when cell sizes are 
                #    are not the same, but regular for each axis.
                canloaddata = True
            # 
            # -- here let's do a simple evaluation on the 
            #    number of cells
            nbins = np.zeros(shape=(3,))
            for ib in range(3):
                if self.binsize[ib] == 0:
                    continue
                # -- follows the x,y,z convention
                nbins[ib] = int(self.domainsize[ib]/self.binsize[ib])
            # 
            # -- get the modflow grid size
            nlay = self.parent.flowmodel.modelgrid.nlay
            nrow = self.parent.flowmodel.modelgrid.nrow
            ncol = self.parent.flowmodel.modelgrid.ncol
            # 
            # -- validate same grid size.
            #    if binsize is zero, then allow
            #    loading if the modflow grid has only
            #    one element in the given dimension
            # 
            # -- x 
            if (nbins[0] != 0):
                if ( nbins[0] != ncol ): 
                    canloaddata = False
            else:
                if ( ncol != 1 ): 
                    canloaddata = False
            # 
            # -- y
            if (nbins[1] != 0):
                if ( nbins[1] != nrow ): 
                    canloaddata = False
            else:
                if ( nrow != 1 ): 
                    canloaddata = False
            # 
            # -- z
            if (nbins[2] != 0):
                if ( nbins[2] != nlay ): 
                    canloaddata = False
            else:
                if ( nlay != 1 ): 
                    canloaddata = False
            #
            # -- load the gpkde data into an array
            if canloaddata:
                # 
                data = np.empty((ntimes, nlay, nrow, ncol), dtype=np.float32)
                data[:,:,:,:] = defaultnan
                #
                # -- load for all times
                for it, time in enumerate(self.times):
                    # -- get data for this time
                    tdata = filtdata[filtdata['time'] == time]
                    #
                    # -- assign spatial array
                    layers = np.unique( tdata['idbinz'] )
                    for lay in layers:
                        srec = tdata[ tdata['idbinz'] == lay ]
                        laydata = np.empty((nrow,ncol),dtype=np.float32)
                        laydata[:,:] = defaultnan
                        laydata[ nrow - srec['idbiny'] - 1, srec['idbinx'] ] = srec[which] 
                        data[it, nlay - lay - 1,:,:] = laydata
                # -- return
                return data
            else:
                print( 
                    f"Warning: data for flowmodel with non-regular StructuredGrid" 
                    f" is returned filtered only by speciesid. The user should apply"
                    f" an adequate coordinates conversion (e.g., scipy.interpolate.griddata)."
                )
                # -- return
                return filtdata
        else:
            print( 
                f"Warning: data for flowmodel with {str(type(self.parent.flowmodel.modelgrid))} grid"
                f" is returned filtered only by speciesid. The user should apply"
                f" an adequate coordinates conversion. "
            )
            # return
            return filtdata

