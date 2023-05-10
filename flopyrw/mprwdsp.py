'''
Configuration of MODPATH-RW dispersion
'''

# python
import numpy as np

# flopy
from flopy.pakbase import Package
from flopy.utils import Util3d

# local 
from .utils import count_instances # Increment COUNTER


class ModpathRWDsp( Package ):
    """
    MODPATH-RW Dispersion Package Class.
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
        modelkind         = None , # linear, nonlinear
        modelkindid       = None , # 1:linear, 2:nonlinear
        alphal            = 0.1  , # xx
        alphath           = 0.01 , # yy
        alphatv           = 0.01 , # zz
        dmeff             = 0.0  , # corrected by tortuosity
        dmaqueous         = 0.0  ,
        betal             = 1    ,
        betath            = 0.5  , 
        betatv            = 0.5  , 
        delta             = 5    , 
        dgrain            = 1    ,
        id                = None , 
        stringid          = None ,
        extension         = 'dsp',
    ):

        # Define UNITNUMBER if the first instance is created
        if self.__class__.COUNTER == 1:
            unitnumber = model.next_unit()
            super().__init__(model, extension, "DSP", unitnumber)
            self.__class__.UNITNUMBER = self.unit_number[0]
        else:
            # Needed for Util3d
            self._parent = self.INSTANCES[0]._parent
        

        # What about this ?
        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])
        self.model = model
        self.shape3d = shape3d


        # Select modelkind
        if modelkind is not None: 
            if modelkind not in ('linear', 'nonlinear'):
                raise Exception( 'flopyrwpt.py: modelkind ' + modelkind + ' is not allowed: linear or nonlinear' )
            if modelkind == 'linear':
                self.modelkind   = modelkind
                self.modelkindid = 1
            elif modelkind == 'nonlinear':
                self.modelkind   = modelkind
                self.modelkindid = 2
        else:
            self.modelkind   = 'linear'
            self.modelkindid = 1


        # Read dispersion model parameters 
        if self.modelkindid == 1:
            # linear
            self.alphal = Util3d(
                model,
                shape3d,
                np.float32,
                alphal,
                name="ALPHAL",
                locat=self.__class__.UNITNUMBER,
            )
            self.alphath = Util3d(
                model,
                shape3d,
                np.float32,
                alphath,
                name="ALPHATH",
                locat=self.__class__.UNITNUMBER,
            )
            self.alphatv = Util3d(
                model,
                shape3d,
                np.float32,
                alphatv,
                name="ALPHATV",
                locat=self.__class__.UNITNUMBER,
            )
            # Consider some checks to avoid writing unnecessary zeroes
            self.dmeff = Util3d(
                model,
                shape3d,
                np.float32,
                dmeff,
                name="DMEFF",
                locat=self.__class__.UNITNUMBER,
            )
        elif self.modelkindid == 2:
            # nonlinear
            self.dmaqueous = dmaqueous
            # Consider some checks to avoid writing unnecessary zeroes
            self.dmeff = Util3d(
                model,
                shape3d,
                np.float32,
                dmeff,
                name="DMEFF",
                locat=self.__class__.UNITNUMBER,
            )
            self.betal = Util3d(
                model,
                shape3d,
                np.float32,
                betal,
                name="BETAL",
                locat=self.__class__.UNITNUMBER,
            )
            self.betath= Util3d(
                model,
                shape3d,
                np.float32,
                betath,
                name="BETATH",
                locat=self.__class__.UNITNUMBER,
            )
            self.betatv= Util3d(
                model,
                shape3d,
                np.float32,
                betatv,
                name="BETATV",
                locat=self.__class__.UNITNUMBER,
            )
            self.delta = Util3d(
                model,
                shape3d,
                np.float32,
                delta,
                name="DELTA",
                locat=self.__class__.UNITNUMBER,
            )
            self.dgrain = Util3d(
                model,
                shape3d,
                np.float32,
                dgrain,
                name="DGRAIN",
                locat=self.__class__.UNITNUMBER,
            )



        # Define id
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER

        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = 'DSP'+str(self.__class__.COUNTER)


        # Add package and save instance        
        if self.__class__.COUNTER == 1:
            self.parent.add_package(self)
        self.__class__.INSTANCES.append( self )


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
            
            # Write id's
            f.write(f"{ins.id}\n")
            f.write(f"{ins.stringid}\n")

            # Write modelkindid
            f.write(f"{ins.modelkindid}\n")
           
            # Write dispersion model parameters
            if ins.modelkindid == 1:
                f.write(ins.dmeff.get_file_entry())
                f.write(ins.alphal.get_file_entry())
                f.write(ins.alphath.get_file_entry())
                f.write(ins.alphatv.get_file_entry())

            if ins.modelkindid == 2:
                f.write(f"{ins.dmaqueous:.16f}\n")
                f.write(ins.dmeff.get_file_entry())
                f.write(ins.betal.get_file_entry())
                f.write(ins.betath.get_file_entry())
                f.write(ins.betatv.get_file_entry())
                f.write(ins.delta.get_file_entry())
                f.write(ins.dgrain.get_file_entry())
       

        # Done
        f.close()

        return
