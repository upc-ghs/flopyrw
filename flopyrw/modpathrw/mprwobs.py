'''
Configuration of MODPATH-RW observations 
'''

# python
import numpy as np
from enum import Enum

# flopy
from flopy.pakbase import Package
from flopy.utils import Util3d
from flopy.discretization import StructuredGrid

# local 
from .utils import multipackage


class observationKindOption(Enum): 
    '''
    Enumerate formats for observation kinds
    '''
    res      = 0
    resident = 0
    flux     = 1
    sink     = 1


class histogramBinOption(Enum): 
    '''
    Enumerate formats for observation kinds
    '''
    scott    = 0
    fd       = 1
    given    = 2
    ts       = 3


class ModpathRWObs( Package ):
    """
    MODPATH-RW Observation Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class: ModpathRW) to which
        this package will be added.
    kind  : int or str
        The kind of observation. Allowed values are:
        * 0 or ('res','resident'): observation of resident concentrations, any cell.
        * 1 or ('flux','sink')   : observation of flux concentrations. For consistency, 
             it should be applied to cells with sink flow and 
             configure the simulation with weaksinkoption=stop_at.
    cellinputoption : int
        The format in which observation cells will be given. Allowed values are:
        * 0: cells are given as a list of individual cells.
        * 1: cells read with free-format input.
    cells : list, np.array
        The list of cells to be considered for the observation.
        * If cellinputoption == 0, then this parameter contain the list of cell ids
        * If cellinputoption == 1, this should an array with the dimensions
          of the flow model. In this last case, array entries containing ones
          are included as part of the observation.
    structured : bool 
        Indicates the format of cell ids in the cells param. If True, then 
        cells are in format (lay, row, col), and if False these are the 
        linear cell numbers.
    outputoption : int
        Defines the information stored in the output file. Allowed values are:
        * 0: output file contains the source particle records
        * 1: postprocessed observation in output file (timeseries) 
             and a second file with source records
        * 2: only one file with the postprocessed observation (timeseries)
    postprocessoption : int
        The kind of postprocess to be applied to the source observation records. 
        Allowed values are:
        * 0: timeseries with histogram reconstruction
        * 1: timeseries with both histogram and smoothed reconstruction
    reconstructionoptions : bool
        Enable interpretation of reconstruction parameters. 
        * False: timeseries reconstruction with default parameters. 
        * True : timeseries reconstruction with parameters provided by the user.
    noptloops : int
        The number of bandwidth optimization loops. It should be a positive integer.
    convergence : float
        The limit relative change of the variables monitored during bandwidth optimization 
        to determine the convergence of reconstruction. It should be positive.
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
        the bandwidth optimization. Amplifies the bin size. In timeseries reconstruction, 
        the bin size is time step defined for the timeseries simulation. 
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
    histogramoptions : bool 
        Enable interpretation of histogram parameters. The histogram is the base over 
        which smoothed reconstruction is performed. 
        * False: histogram options handled internally with default parameters, not written. 
        * True : interpret histogram options, written to package file.
    binoption : int or str
        Selection of the histogram bin size. Allowed values are: 
          * 0 or 'scott': applies the automatic selection for Gaussian distributions.
          * 1 or 'fd'   : applies the automatic selection based on the inter quartile range.
          * 2 or 'given': user provided value in binparam.
          * 3 or 'ts'   : infered from the timeseries, assumed uniform time step.
    binparam : float
        Auxiliary parameter controlling the selection of the histogram bin size. 
        If binoption is 0 or 1, binparam is fraction (0,1] multiplying the automatic value. 
        If binoption is 2, binparam is the bin size. 
        If binoption is 3, binparam is not read.
    timestepout : float
        Value determining the time discretization to write the timeseries data. If given, 
        interpolates the values obtained from the model to the given resolution. 
        Interpolation is linear and performed only within the limits of the real data. 
    filename : str
        The output file where the observation is written, the name 
        without the extension. If None is given, will employ stringid. 
        * If outputoption == 0, this file contains the observation records
        * If outputoption == 1, this file contains the postprocessed 
          observation and a second file whose name is rec+filename
          contains the observations records. 
        * If outputoption == 2, this file contains the postprocessed 
          observation timeseries
    stringid : str
        An identifier for the observation. By default assign the concatenation
        of package ftype and an integer id. 
    extension : str
        The extension for the observation file. By default 'obs' 
    """


    @staticmethod
    def _ftype():
        return 'OBS'


    @multipackage
    def __init__(
        self,
        model,
        kind                   = 0,
        cellinputoption        = 0,
        cells                  = None,
        structured             = None, 
        outputoption           = 2,
        postprocessoption      = 1,
        reconstructionoptions  = False, 
        noptloops              = 10, 
        convergence            = 0.01,
        initialsmoothingformat = 0,
        binsizefactor          = 3.0,
        effectiveweightformat  = 0,
        histogramoptions       = False, 
        binoption              = 1, 
        binparam               = None,
        timestepout            = None, 
        filename               = None, 
        stringid               = None,
        extension              = 'obs',
    ):

        # Call parent constructor
        ftype = self._ftype()
        if ( model.multipackage[ftype]['count'] == 0 ): 
            super().__init__(
                model,
                extension,
                ftype,
                model.multipackage[ftype]['unitnumber']
            )
        else:
            # Pass the parent
            self._parent = model.multipackage[ftype]['instances'][0]._parent


        # Observation kind
        if ( not isinstance( kind, (int,str) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for observation kind."
                f" It can be specified as int or str, but {str(type(kind))}"
                f" was given." 
            )
        if ( isinstance( kind, int ) ):
            if ( kind not in [0,1] ): 
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid kind option {str(kind)}."
                    f" Allowed values are 0 (resident) or 1 (flux)."
                )
            self.kind = kind
            if self.kind == 0: 
                self.stringkind = 'RESIDENT'
            elif self.kind == 1:
                self.stringkind = 'FLUX'
        elif ( isinstance( kind, str ) ):
            try:
                self.kind = observationKindOption[kind.lower()].value
                if self.kind == 0: 
                    self.stringkind = 'RESIDENT'
                elif self.kind == 1:
                    self.stringkind = 'FLUX'
            except:
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid observation kind {str(kind)}."
                    f" The allowed values are resident or flux."
                )

        # output option
        if ( outputoption not in [0,1,2] ):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid output option {str(kind)}."
                f" Allowed values are 0 (only obs records),"
                f" 1 (records and postprocess) or 2 (only postprocess)."
            )
        self.outputoption = outputoption

        # postprocess option
        if ( postprocessoption not in [0,1] ):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid postprocess option {str(kind)}."
                f" Allowed values are 0 (only histogram) or 1"
                f" (histogram and smoothed reconstruction)."
            )
        self.postprocessoption = postprocessoption

        # cellinputoption 
        if ( cellinputoption not in [0,1] ):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid cellinputoption option {str(cellinputoption)}"
                f" Allowed values are 0 (list of cellids) or 1 (modelgrid array)."
            )
        self.cellinputoption = cellinputoption

        if self.cellinputoption == 0:

            if ( cells is None ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid cells specification. It was expecting a list"
                    f" of cell ids, but None was given."
                )

            # Infer structured parameter from 
            # the modelgrid type if it remained with 
            # default None
            if ( structured is None ):
                structured = False
                if ( isinstance( self._parent.flowmodel.modelgrid, StructuredGrid ) ):
                    # If grid is structured, infers True
                    structured = True

            if ( not isinstance( structured, bool ) ): 
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for structured parameter."
                    f" It should be boolean, but"
                    f" {str(type(structured))} was given."
                )

            # Cell format
            if structured:
                self.cellformat = 0
            else:
                self.cellformat = 1
            self.structured = structured
            
            # Some validation with the modelgrid
            if (
                (self.structured) and
                not isinstance( self._parent.flowmodel.modelgrid, StructuredGrid )
            ):
                raise Exception(
                    f"{self.__class__.__name__}:"
                    f" Invalid cellinputoption and structured parameters"
                    f" for the modelgrid of type {str(type(self._parent.flowmodel.modelgrid))}."
                    f" Needs to be StructuredGrid in order to specify cells as (lay,row,col)."
                )

            # Write as a list of cellids
            if cells is None:
                self.cells = []
            else:
                if not isinstance( cells, (list,np.ndarray,tuple) ):
                    raise TypeError(
                        f"{self.__class__.__name__}:"
                        f" Invalid input cells type."
                        f" It should be list, tuple of np.ndarray."
                        f" {str(type(cells))} was given."
                    )
                # Transform to numpy array
                cells = np.array( cells )
                
                if self.structured: 
                    # The case of a single list/tuple with cell id
                    if ( (len(cells.shape)==1) ): 
                        if ( len(cells)==3 ):
                            cells = np.expand_dims( cells, axis=0 )
                        else:
                            raise ValueError(
                                f"{self.__class__.__name__}:"
                                f" Invalid cells specification. While using the"
                                f" structured flag for cellinputoption = 0, the cells"
                                f" should be specified as a list of cell ids following"
                                f" the format (lay,row,col)."
                            )
                    elif ( len(cells.shape)==2 ):
                        if ( not ( cells.shape[1] == 3 ) ): 
                            raise ValueError(
                                f"{self.__class__.__name__}:"
                                f" Invalid cells specification. While using the"
                                f" structured flag for cellinputoption = 0, the cells"
                                f" should be specified as a list of cell ids following"
                                f" the format (lay,row,col)."
                            )
                    else:
                        raise ValueError(
                            f"{self.__class__.__name__}:"
                            f" Invalid cells specification. While using the"
                            f" structured flag for cellinputoption = 0, the cells"
                            f" should be specified as a list of cell ids following"
                            f" the format (lay,row,col)."
                        )
                else:
                    # For unstructured, ask for a one-dimensional array
                    if ( (len(cells.shape)!=1) ): 
                        raise ValueError(
                            f"{self.__class__.__name__}:"
                            f" Invalid cells specification. While using"
                            f" the unstructured format for cellinputoption=0,"
                            f" cell ids should be given as linear indexes format."
                        )
                # Validate dtype
                if ( cells.dtype not in [np.int32,np.int64,int] ):
                    raise TypeError(
                        f"{self.__class__.__name__}:"
                        f" Invalid type for cells specification." 
                        f" Cell ids should be of integer type, but"
                        f" {str(cells.dtype)} was given."
                    )

                # To the class
                self.cells = cells

            # If no cells, stop. 
            if len(self.cells) == 0:
                raise ValueError( 
                    f"{self.__class__.__name__}:"
                    f" No cells were given for the observation."
                )

        elif self.cellinputoption == 1: 
            # Write obs cells as 3D array
            shape = model.shape
            if len(shape) == 3:
                shape3d = shape
            elif len(shape) == 2:
                shape3d = (shape[0], 1, shape[1])
            else:
                shape3d = (1, 1, shape[0])
            self.model   = model
            self.shape3d = shape3d

            if ( cells is None ):
                raise ValueError( 
                    f"{self.__class__.__name__}:"
                    f" No cells were given for the observation."
                )
            else:
                if not isinstance( cells, (int,list,np.ndarray) ):
                    raise TypeError(
                        f"{self.__class__.__name__}:"
                        f" Invalid type for input cells. It should be list or np.ndarray."
                        f" {str(type(cells))} was given."
                    )
                cells  = np.array(cells)
                ucells = np.unique(cells)
                for uc in ucells: 
                    if uc not in [0,1]:
                        raise ValueError( 
                            f"{self.__class__.__name__}:"
                            f" Invalid values for cells specification."
                            f" While using cellinputoption=1, the cells specification"
                            f" should only contain values 0 or 1."
                            f" The value {str(uc)} was found."
                        )

                try:
                    self.cells = Util3d(
                        model,
                        shape3d,
                        np.int32,
                        cells,
                        name="OBSCELLS",
                        locat=self._parent.multipackage[ftype]['unitnumber'], 
                    )
                except:
                    raise Exception(
                        f"{self.__class__.__name__}:"
                        f" Error while initializing distributed variable cells."
                        f" Is the input shape consistent with flow model dimensions ?"
                    )

        if ( not isinstance( reconstructionoptions,bool ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for reconstructionoptions. It should be bool, but"
                f" {str(type(reconstructionoptions))} was given."
            )
        if reconstructionoptions:
            self.reconstructionoptions = 1
        else:
            self.reconstructionoptions = 0

        # process reconstruction options
        if reconstructionoptions:

            # noptloops
            if ( not isinstance( noptloops,int ) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for noptloops. It should be int, but"
                    f" {str(type(noptloops))} was given."
                )
            if ( noptloops < 0 ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid value for noptloops. It should be positive or zero."
                )
            self.noptloops = noptloops

            # convergence
            if ( not isinstance( convergence, (int,float) ) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for convergence. It should be float, but"
                    f" {str(type(convergence))} was given."
                )
            if ( convergence <= 0 ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid value for convergence. It should be a positive real number."
                )
            self.convergence = convergence

            # initialsmoothingformat
            if ( not isinstance( initialsmoothingformat,int ) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for initialsmoothingformat. It should be int, but"
                    f" {str(type(initialsmoothingformat))} was given."
                )
            if initialsmoothingformat not in [0,1]: 
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid value for initialsmoothinformat. Allowed values are"
                    f" 0 (automatic for Gaussian) or 1 (factor multiplying bin size)."
                )
            self.initialsmoothingformat = initialsmoothingformat

            # binsizefactor
            if ( not isinstance( binsizefactor, (int,float) ) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for binsizefactor. It should be float, but"
                    f" {str(type(binsizefactor))} was given."
                )
            if ( binsizefactor <= 0 ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid value for binsizefactor. It should be a positive real number."
                )
            self.binsizefactor = binsizefactor

            # effectiveweightformat
            if ( not isinstance( effectiveweightformat,int ) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for effectiveweightformat. It should be int, but"
                    f" {str(type(effectiveweightformat))} was given."
                )
            if effectiveweightformat not in [0,1,2,3]: 
                raise ValueError(
                    f"{self.__class__.__name__}:" 
                    f" Invalid value for effectiveweightformat. Determines transformation"
                    f" to an equivalen count of particles. Allowed values are"
                    f" 0 (unique effective mass), 1 (scale by averaged mass),"
                    f" 2 (real histogram) or 3 (local effective histogram)."
                    f" {str(effectiveweightformat)} was given."
                )
            self.effectiveweightformat = effectiveweightformat


        if ( not isinstance( histogramoptions,bool ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for histogramoptions. It should be bool, but"
                f" {str(type(histogramoptions))} was given."
            )
        if histogramoptions:
            self.histogramoptions = 1
        else:
            self.histogramoptions = 0

        if histogramoptions: 

            # binoption
            if ( not isinstance( binoption,int ) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for binoption. It should be int, but"
                    f" {str(type(binoption))} was given."
                )
            if binoption not in [0,1,2,3]: 
                raise ValueError(
                    f"{self.__class__.__name__}:" 
                    f" Invalid value for binoption. Determines selection"
                    f" of histogram bin size. Allowed values are"
                    f" 0 (based on Scott's rule), 1 (based on Freedman-Diaconis rule),"
                    f" 2 (given by user) or 3 (infered from timepointdata)."
                    f" {str(binoption)} was given."
                )
            self.binoption = binoption

            # binparam
            if binparam is None: 

                if ( self.binoption in [0,1] ):
                    # Set a default 
                    self.binparam = 0.75

                if ( self.binoption == 2 ):
                    # Needs to be requested
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Invalid value for binparam. None was given."
                    )

            else:

                if ( not isinstance( binparam, (int,float) ) ):
                    raise TypeError(
                        f"{self.__class__.__name__}:"
                        f" Invalid type for binparam. It should be float, but"
                        f" {str(type(binparam))} was given."
                    )

                if ( self.binoption in [0,1] ):
                    if ( ( binparam <= 0 ) or ( binparam > 1 ) ):
                        raise ValueError(
                            f"{self.__class__.__name__}:"
                            f" Invalid value for binparam. It should be a fraction."
                        )
                    self.binparam = binparam

                if ( self.binoption == 2 ):
                    if ( binparam <= 0 ):
                        raise ValueError(
                            f"{self.__class__.__name__}:"
                            f" Invalid value for binparam. It should be a positive real number."
                        )
                    self.binparam = binparam

            # timestepout
            if timestepout is not None: 
                if ( not isinstance( timestepout, (int,float) ) ):
                    raise TypeError(
                        f"{self.__class__.__name__}:"
                        f" Invalid type for timestepout. It should be float, but"
                        f" {str(type(timestepout))} was given."
                    )
                if ( timestepout <= 0 ):
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Invalid value for timestepout. It should be a positive real number."
                    )
                self.timestepout = timestepout
            else:
                self.timestepout = None


        # String id
        if (stringid is not None):
            if ( not isinstance( stringid, str ) ): 
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for stringid. It should be str, but"
                    f" {str(type(stringid))} was given."
                )
            self.stringid = stringid
        else:
            self.stringid = ftype+str(self.id) 

        # filename   
        if ( filename is None ):
            # Assign by default the string id, in lower case
            filename = self.stringid.lower()
        else:
            # If a value was given, validate
            if ( not isinstance( filename, str ) ): 
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for filename. It should be str, but"
                    f" {str(type(filename))} was given."
                )
        # Some validation for extension
        if ( not isinstance( extension, str ) ): 
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for extension. It should be str, but"
                f" {str(type(extension))} was given."
            )
        # Filename for this observation
        self.filename = f"{filename}.{extension}"

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

        # Shortcut to multipackage dict
        mpkg = self.parent.multipackage[self._ftype()]
        
        # Write
        with open( mpkg['instances'][0].fn_path, 'w' ) as f:

            # Write how many instances...
            f.write(f"{mpkg['count']}\n")

            # ...and loop over them
            for ins in mpkg['instances']:

                # Write obs stringid
                f.write(f"{ins.stringid}\n")

                # Write obs stringkind
                f.write(f"{ins.stringkind}\n")

                # Write obs filename
                f.write(f"{ins.filename}\n")

                # For sink observations
                if ins.kind == 1:
                    # Write the observation options 
                    f.write(f"{ins.outputoption}  {ins.postprocessoption}  {ins.histogramoptions}  {ins.reconstructionoptions}\n")
                else:
                    # Write the observation options 
                    f.write(f"{ins.outputoption}  {ins.postprocessoption}  {ins.reconstructionoptions}\n")

                # write the histogram data
                if ( (ins.kind==1) and (ins.histogramoptions==1) ):
                    if ( ins.binoption==3 ): 
                        f.write(f"{ins.binoption}\n")
                    else:
                        if ( ins.timestepout is None ): 
                            f.write(f"{ins.binoption}  {ins.binparam:.10f}\n")
                        else:
                            f.write(f"{ins.binoption}  {ins.binparam:.10f}  {ins.timestepout:.10f}\n")

                if ins.reconstructionoptions == 1:
                    # Write noptloops and convergence
                    f.write(f"{ins.noptloops}   {ins.convergence}\n")

                    if ins.initialsmoothingformat == 1:
                        # Write initialsmoothingformat and binsizefactor
                        f.write(f"{ins.initialsmoothingformat}   {ins.binsizefactor}\n")
                    else:
                        # Write initialsmoothingformat
                        f.write(f"{ins.initialsmoothingformat}\n")

                    # Write effectiveweightformat
                    f.write(f"{ins.effectiveweightformat}\n")



                # Write the celloption params
                if ins.cellinputoption == 0:
                    # The flag and specs
                    f.write(f"{ins.cellinputoption}  {len(ins.cells)}  {ins.cellformat} \n")

                    # And the list
                    fmts = []
                    if ins.structured:
                        fmts.append("{:9d}") # lay
                        fmts.append("{:9d}") # row
                        fmts.append("{:9d}") # col
                    else:
                        fmts.append("{:9d}") # cellid
                    fmt = " " + " ".join(fmts) + "\n"
                    for oc in ins.cells:
                        woc = np.array(oc).astype(np.int32)+1 # Correct the zero-based indexes
                        if ins.structured:
                            f.write(fmt.format(*woc))
                        else:
                            f.write(fmt.format(woc))

                elif ins.cellinputoption == 1:

                    # The flag
                    f.write(f"{ins.cellinputoption}\n")

                    # And the array
                    f.write(ins.cells.get_file_entry())


        # Done
        return


    def get_output(self):
        '''
        Load the observation output file into a rec array

        Note: the current implementation loads data for the cases
              where the time vector is different for each observation, 
              that is, the bin size was selected automatically and without
              interpolation (data appended vertically).
        '''
        import os 
        from flopy.utils.flopy_io import loadtxt
        # Shouldn't this be a case to handle at loadtxt ?
        try:
            import pandas as pd
            use_pandas = True
        except ImportError:
            use_pandas = False
        

        if self.stringkind == 'RESIDENT':
            raise NotImplementedError(
                    f"{self.__class__.__name__}:get_output:"
                    f" Output loading is not yet implemented for observations of kind RESIDENT."
                )

        if self.stringkind == 'FLUX':
            if self.outputoption in [1,2]:
                # Loading structure at this point only for btc's
                # with specific time vector
                # This needs some review: in simulations 
                # with adaptgridtocoords = true (an endpoint for example)
                # the time vector can also be different for each species, 
                # although interpolation was applied.
                if (
                    ( self.histogramoptions==1 ) and (
                        ( self.binoption in [2,3] ) 
                        #or ( self.timestepout is not None )
                    )
                ):
                    raise NotImplementedError(
                            f"{self.__class__.__name__}:get_output:"
                            f" Output loading is not yet implemented for observations of type FLUX "
                            f" with the same time step for all species."
                        )

                # continue
                dtypelist = [
                        ("speciesid" , np.int32  ),
                        ("tid"       , np.int32  ),
                        ("time"      , np.float32),
                        ("qsink"     , np.float32),
                        ("chist"     , np.float32),
                    ]
                # extend list if postprocessed
                if self.postprocessoption == 1:
                    dtypelist.append( ("cgpkde"    , np.float32) )

                # load
                recdata = loadtxt(
                    os.path.join( self._parent.model_ws, self.filename ),
                    dtype=np.dtype(dtypelist),
                    skiprows=0,
                    use_pandas=use_pandas
                )
                
                # To python zero-based indexation
                recdata['speciesid'] = recdata['speciesid'] - 1

                # done
                return recdata

            # for observations with only records
            else:
                raise NotImplementedError(
                        f"{self.__class__.__name__}:get_output:"
                        f" Output loading is not yet implemented for observations of type FLUX "
                        f" storing only records."
                    )
