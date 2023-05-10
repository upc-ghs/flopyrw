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
        kind              = 0,
        outputoption      = 2,
        postprocessoption = 1,
        cellinputoption   = 0,
        cells             = None,
        structured        = True, 
        basefilename      = 'mpathrwobs_',
        filename          = None,
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
