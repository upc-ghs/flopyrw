'''
Configuration of MODPATH-RW initial conditions
'''

# python
import numpy as np

# flopy
from flopy.pakbase import Package
from flopy.utils import Util3d

# local 
from .utils import count_instances # Increment COUNTER


class ModpathRWIc( Package ):
    """
    MODPATH-RW Initial Condition Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class: ModpathRW) to which
        this package will be added.
    kind  : int 
        The kind of initial condition specification. In the meantime, there is only 
        one format enabled loading the concentration parameter with u3d.
    mass  : float
        The scale of particles mass used for transforming the concentration data into 
        a distribution of particles. 
    concentration : list, ndarray, int or float
        Interpreted as a resident initial concentration. It will be written with Util3d.
    speciesid : int 
        The species id to which this initial concentration is related. Will only 
        be written if particlesmassoption == 2. It should be given as the zero-based 
        index of the specie (the python-index of the specie in the list of species), 
        but is written to the package file with fortran format (+1, see the write_file method).
    id : int 
        A positive integer identifier. If None given is automatically assigned with 
        the instance counter. 
    stringid : str
        A string identifier. If None is given is automatically generated as IC+str(id).
    extension : str
        The file extension for the package. Defaults to ic. 
    """


    # Class properties
    COUNTER    = 0
    UNITNUMBER = 0
    INSTANCES  = []


    @count_instances
    def __init__(
        self,
        model,
        kind          = 0   ,
        mass          = 1.0 ,
        concentration = None,
        speciesid     = None,
        id            = None,
        stringid      = None,
        extension     = 'ic',
    ):

        # Define UNITNUMBER if the first instance is created
        if self.__class__.COUNTER == 1:
            unitnumber = model.next_unit()
            super().__init__(model, extension, "IC", unitnumber)
            self.__class__.UNITNUMBER = self.unit_number[0]
        else:
            # Needed for Util3d
            self._parent = self.INSTANCES[0]._parent

        # Determines format for writing the initial condition
        if (kind not in [0]): 
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid kind parameter ' +str(kind)+
                '. Allowed values are 0 (read concentration with u3d).'
            )
        self.kind = kind

        # The particles mass
        if (mass <= 0):
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid particle mass for initial condition. ' +
                'It should be a positive value.' 
            )
        self.mass = mass

        # A species id, will only be written 
        # if particlesmassoption == 2
        self.speciesid = speciesid

        # define id
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER

        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = 'IC'+str(self.id)

        # and the concentration: needs id
        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])
        self.model = model
        self.shape3d = shape3d

        if (
            ( concentration is not None ) and 
            ( self.kind == 0 )
        ):
            # TypeError if any other unsupported type
            if ( not isinstance( concentration, (list,float,int,np.ndarray) ) ):
                raise TypeError(
                    self.__class__.__name__ + ':' + 
                    ' Invalid type ' + str(type(concentration)) + 
                    ' for concentration. It should be list, np.ndarray or a single numeric value.'
                )
            # If list, to np.array
            if ( isinstance( concentration, list ) ): 
                concentration = np.array( concentration )
            # Verify a valid dtype
            if ( isinstance( concentration, np.ndarray ) ): 
                if ( 
                    (concentration.dtype not in [np.int32, np.int64, np.float32, np.float64])
                ):
                    raise TypeError(
                        self.__class__.__name__ + ':' + 
                        ' Invalid dtype ' + str(concentration.dtype) + 
                        ' for concentration. It should be int or float. ' 
                    )
            self.concentration = Util3d(
                model,
                shape3d,
                np.float32,
                concentration,
                name=self.stringid,
                locat=self.__class__.UNITNUMBER,
            )

        elif ( (self.kind == 0) and (concentration is None) ):
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid initial concentration distribution. None was given.'
            )
     

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

            # Write string id 
            f.write(f"{ins.stringid}\n")

            # Kind/format of initial condition 
            f.write(f"{ins.kind}\n")

            # 1: resident concentration array  
            if ins.kind == 0:

                # Give particles mass 
                f.write(f"{ins.mass:.10f}\n")

                # If solute ids are being specified, write the id
                if( 
                    ( self.INSTANCES[0]._parent.particlesmassoption == 2 ) 
                ):
                    if ( 
                        ( ins.speciesid is None ) or
                        ( ins.speciesid < 0 ) 
                    ):
                        raise ValueError(
                            self.__class__.__name__ + ':' + 
                            ' Invalid species id for initial condition. ' +
                            'It should be a valid index in the list of species.' 
                        )
                    f.write(f"{ins.speciesid+1}\n")


                # And the concentration distribution
                f.write(ins.concentration.get_file_entry())


        # Done
        return
