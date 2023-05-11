'''
Configuration of MODPATH-RW observations 
'''

# python
import numpy as np

# flopy
from flopy.pakbase import Package
from flopy.utils import Util3d

# local 
from .utils import count_instances # Increment COUNTER


class ModpathRWObs( Package ):
    """
    MODPATH-RW Observation Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class: ModpathRW) to which
        this package will be added.
    kind  : int
        The kind of observation. Allowed values are:
            0: observation of resident concentrations, any cell
            1: observation of flux concentrations. For consistency, it should be applied 
               to cells with sink flow and with weaksinkoption is stop_at
    cellinputoption : int
        The format in which observation cells will be given. Allowed values are:
            0: cells are given as a list of individual cells
            1: cells read with free-format input
    cells : list, np.array
        The list of cells to be considered for the observation.
        If cellinputoption == 0, then this parameter contain the list of cell ids
        If cellinputoption == 1, this should an array with the dimensions of the flow model. 
        In this last case, array entries containing ones are included as part of the observation.
    structured : bool 
        Indicates the format of cell ids in the cells param. If True, then cells are in format (lay, row, col), 
        and if False these are the linear cell numbers.
    outputoption : int
        Defines the information stored in the output file. Allowed values are:
            0: output file contains the source particle records
            1: postprocessed observation in output file (timeseries) and a second file with source records
            2: only one file with the postprocessed observation (timeseries)
    postprocessoption : int
        The kind of postprocess to be applied to the source observation records. Allowed values are:
            0: timeseries with histogram reconstruction
            1: timeseries with both histogram and smoothed reconstruction
    basefilename : str
        The output file where the observation is written. 
        If outputoption == 0, this file contains the observation records
        If outputoption == 1, this file contains the postprocessed observation and a second file
        whose name is rec+basefilename contains the observations records. 
        If outputoption == 2, this file contains the postprocessed observation timeseries
    id : int 
        A positive integer identifier. If None given is automatically assigned with 
        the instance counter. 
    stringid : str
        An identifier for the observation. By default the name 'OBS'+str(id) is assigned
    extension : str
        The extension for the observation file. By default 'obs' 
    """

    # Class properties
    COUNTER    = 0
    UNITNUMBER = 0
    INSTANCES  = []
 
    @count_instances
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
        id                = None,
        stringid          = None,
        extension         = 'obs',
    ):

        # Define UNITNUMBER for the first instance
        if self.__class__.COUNTER == 1:
            unitnumber = model.next_unit()
            super().__init__(model, extension, "OBS", unitnumber)
            self.__class__.UNITNUMBER = self.unit_number[0]
        else:
            # Needed for Util3d
            self._parent = self.INSTANCES[0]._parent

        # Observation kind
        if ( kind not in [0,1] ): 
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid kind option ' +str(kind)+
                '. Allowed values are 0 (Resident) or 1 (Flux).'
            )
        self.kind = kind
        if self.kind == 0: 
            self.stringkind = 'RESIDENT'
        elif self.kind == 1:
            self.stringkind = 'FLUX'

        # output option
        if ( outputoption not in [0,1,2] ):
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid output option ' +str(kind)+
                '. Allowed values are 0 (only obs records), 1 (records and postprocess) or 2 (only postprocess).'
            )
        self.outputoption = outputoption

        # postprocess option
        if ( postprocessoption not in [0,1] ):
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid postprocess option ' +str(kind)+
                '. Allowed values are 0 (only histogram) or 1 (histogram and smoothed reconstruction).'
            )
        self.postprocessoption = postprocessoption

        # cellinputoption 
        if ( cellinputoption not in [0,1] ):
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid cellinputoption option ' +str(cellinputoption)+
                '. Allowed values are 0 (list of cellids) or 1 (modelgrid array).'
            )
        self.cellinputoption = cellinputoption

        if self.cellinputoption == 0:
            # Write as a list of cellids
            if cells is None:
                self.cells = []
            else:
                if not isinstance( cells, (list,np.ndarray,tuple) ):
                    raise TypeError(
                        self.__class__.__name__ + ':' + 
                        ' Invalid cells type. It should be list, tuple of np.array. ' + 
                        str(type( cells )) + ' was given.'
                    )
                if ( isinstance( cells, tuple ) and len(cells)==1 ):
                    cells = cells[0]
                # Maybe some sanity check about data structure or the same 
                # used for partlocs
                self.cells = cells

            # If no cells, stop. 
            if len(self.cells) == 0:
                raise ValueError( 
                    self.__class__.__name__ + ':' + 
                    ' No cells were given for the observation.'
                )

            # Cell format
            if structured:
                self.cellformat = 0
            else:
                self.cellformat = 1
            self.structured = structured

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
                    self.__class__.__name__ + ':' + 
                    ' No cells were given for the observation.'
                )
            else:
                if not isinstance( cells, (list,np.ndarray) ):
                    raise TypeError(
                        self.__class__.__name__ + ':' + 
                        ' Invalid cells type. It should be list, or np.array. ' + 
                        str(type( cells )) + ' was given.'
                    )
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
