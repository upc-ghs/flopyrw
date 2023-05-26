'''
Configuration of MODPATH-RW observations 
'''

# python
import numpy as np
from enum import Enum

# flopy
from flopy.pakbase import Package
from flopy.utils import Util3d

# local 
from .utils import multipackage


class observationKindOption(Enum): 
    '''
    Enumerate formats for timeseries output option
    '''
    res      = 0
    resident = 0
    flux     = 1
    sink     = 1


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
    basefilename : str
        The output file where the observation is written. 
        * If outputoption == 0, this file contains the observation records
        * If outputoption == 1, this file contains the postprocessed 
          observation and a second file whose name is rec+basefilename
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
        kind              = 0,
        cellinputoption   = 0,
        cells             = None,
        structured        = True, 
        outputoption      = 2,
        postprocessoption = 1,
        basefilename      = 'mprwobs_',
        stringid          = None,
        extension         = 'obs',
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
                    self.__class__.__name__ + ':'
                    + ' Invalid kind option ' + str(kind)
                    + '. Allowed values are 0 (resident) or 1 (flux).'
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
                self.__class__.__name__ + ':' 
                + ' Invalid output option ' + str(kind)
                + '. Allowed values are 0 (only obs records),'
                + ' 1 (records and postprocess) or 2 (only postprocess).'
            )
        self.outputoption = outputoption

        # postprocess option
        if ( postprocessoption not in [0,1] ):
            raise ValueError(
                self.__class__.__name__ + ':'
                + ' Invalid postprocess option ' + str(kind)
                + '. Allowed values are 0 (only histogram) or '
                + '1 (histogram and smoothed reconstruction).'
            )
        self.postprocessoption = postprocessoption

        # cellinputoption 
        if ( cellinputoption not in [0,1] ):
            raise ValueError(
                self.__class__.__name__ + ':'
                + ' Invalid cellinputoption option '
                + str(cellinputoption)
                + '. Allowed values are 0 (list of cellids) or 1 (modelgrid array).'
            )
        self.cellinputoption = cellinputoption

        if self.cellinputoption == 0:

            if ( not isinstance( structured, bool ) ): 
                raise TypeError(
                    self.__class__.__name__ + ':'
                    + ' Invalid type for structured parameter. '
                    + 'It should be boolean, but ' 
                    + str(type(structured)) + ' was given.'
                )

            # Cell format
            if structured:
                self.cellformat = 0
            else:
                self.cellformat = 1
            self.structured = structured

            # Write as a list of cellids
            if cells is None:
                self.cells = []
            else:
                if not isinstance( cells, (list,np.ndarray,tuple) ):
                    raise TypeError(
                        self.__class__.__name__ + ':'
                        + ' Invalid cells type. '
                        + 'It should be list, tuple of np.ndarray. '
                        + str(type(cells)) + ' was given.'
                    )
                # Transform to numpy array
                cells = np.array( cells )

                if self.structured: 
                    # The case of a single list/tuple with cell id
                    if ( (len(cells.shape)==1) ): 
                        if ( len(cells)==3 ):
                            cells = cells.transpose()
                        else:
                            raise ValueError(
                                self.__class__.__name__ + ':' 
                                + ' Invalid cells specification. While using the '
                                + 'structured flag for cellinputoption = 0, the cells ' 
                                + 'should be specified as a list of cell ids following '
                                + 'the format (lay,row,col).'
                            )
                    elif ( len(cells.shape)==2 ):
                        if ( not ( cells.shape[1] == 3 ) ): 
                            raise ValueError(
                                self.__class__.__name__ + ':' 
                                + ' Invalid cells specification. While using the '
                                + 'structured flag for cellinputoption = 0, the cells ' 
                                + 'should be specified as a list of cell ids following '
                                + 'the format (lay,row,col).'
                            )
                    else:
                        raise ValueError(
                            self.__class__.__name__ + ':' 
                            + ' Invalid cells specification. While using the '
                            + 'structured flag for cellinputoption = 0, the cells ' 
                            + 'should be specified as a list of cell ids following '
                            + 'the format (lay,row,col).'
                        )
                else:
                    # For unstructured, ask for a one-dimensional array
                    if ( (len(cells.shape)!=1) ): 
                        raise ValueError(
                            self.__class__.__name__ + ':' 
                            + ' Invalid cells specification. While using '
                            + 'the unstructured format for cellinputoption=0, ' 
                            + 'cell ids should be given as linear indexes format.'
                        )
                # Validate dtype
                if ( cells.dtype not in [np.int32,np.int64,int] ):
                    raise TypeError(
                        self.__class__.__name__ + ':' 
                        + ' Invalid type for cells specification. ' 
                        + 'Cell ids should be of integer type, but '
                        + str(cells.dtype) + ' was given.'
                    )

                # To the class
                self.cells = cells

            # If no cells, stop. 
            if len(self.cells) == 0:
                raise ValueError( 
                    self.__class__.__name__ + ':'
                    + ' No cells were given for the observation.'
                )

        elif self.cellinputoption == 1: 
            # Write obs cells as 3D array

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
                raise ValueError( 
                    self.__class__.__name__ + ':'
                    + ' No cells were given for the observation.'
                )
            else:
                if not isinstance( cells, (int,list,np.ndarray) ):
                    raise TypeError(
                        self.__class__.__name__ + ':'
                        + ' Invalid type for cells . It should be list, or np.ndarray. '
                        + str(type(cells)) + ' was given.'
                    )
                cells  = np.array(cells)
                ucells = np.unique(cells)
                for uc in ucells: 
                    if uc not in [0,1]:
                        raise ValueError( 
                            self.__class__.__name__ + ':'
                            + ' Invalid values for cells specification. '
                            + 'While using cellinputoption=1, the cells specification'
                            + 'should only contain values 0 or 1. '
                            + 'The value ' + str(uc) + ' was found.'
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
                        self.__class__.__name__ + ':'
                        + ' Error while initializing distributed variable cells. '
                        + 'Is the input shape consistent with flow model dimensions ? '
                    )

        # String id
        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = ftype+str(self.id) 

        # Some validation for basefilename
        if ( not isinstance( basefilename, str ) ): 
            raise TypeError(
                self.__class__.__name__ + ':'
                + ' Invalid type for basefilename. It should be str, but ' 
                + str(type(basefilename)) + ' was given.'
            )
        # Some validation for extension
        if ( not isinstance( extension, str ) ): 
            raise TypeError(
                self.__class__.__name__ + ':'
                + ' Invalid type for extension. It should be str, but ' 
                + str(type(extension)) + ' was given.'
            )
        # Filename for this observation
        self.filename = basefilename+str(self.id)+'.'+extension


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
                
                # Write the output and postprocess options
                f.write(f"{ins.outputoption}  {ins.postprocessoption}\n")

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
