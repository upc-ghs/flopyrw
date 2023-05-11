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
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    binsize : list [sx,sy,sz]
        List with cell size employed for reconstruction process 
    domainsize : list [sx,sy,sz]
        List with domain sizes. These are used for computing reconstruction 
        grid dimensions using given binsize.
    domainorigin : list [ox,oy,oz]
        List with domain origin. 
    noptloops: int
        The number of optimization loops smoothed reconstruction
    asconcentration: bool
        Flag to indicate whether reconstruction should be returned as resident concentration. 
        Because the reconstruction grid is different than flowmodel grid, said transformation 
        can be easily achieved in case both porosities and retardation are spatially uniform.
        If not, then some kind of intersection or relation should be established between the 
        base flowmodel and the reconstruction grid.
    outputfilename: str
        Filename that will be used for writing reconstruction output. Contains 
        cell id, density and both time step and particle group identification.
    extension : str, optional
        File extension (default is 'dispersion').
    """

    def __init__(
        self,
        model,
        binsize         = [1,1,1],
        domainsize      = [1,1,1],
        domainorigin    = None, 
        minhlambda      = 1.0 ,
        maxhlambda      = 0.1 ,
        deltahlambda    = 10.0,
        kerneldatabase  = False, 
        noptloops       = 2,
        asconcentration = False,
        outputfilename  = None,
        extension       = 'gpkde',
    ):

        unitnumber = model.next_unit()
        super().__init__(model, extension, "GPKDE", unitnumber)

        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])

        # Some health checks
        self.binsize    = np.array(binsize).astype(np.float32)
        self.domainsize = np.array(domainsize).astype(np.float32)

        # Fills domain origin with some values 
        # detemined from the modelgrid definition
        if domainorigin is None :
            domainorigin = [
                    self._parent.flowmodel.modelgrid.xoffset,
                    self._parent.flowmodel.modelgrid.yoffset,
                    np.min(self._parent.flowmodel.modelgrid.botm)
                ]
        self.domainorigin = np.array(domainorigin).astype(np.float32)


        # Depending on the kind of modelgrid maybe some
        # decisions can be made for binsizes and/or domain sizes
        # StructuredGrid has the is_regular function

        self.noptloops       = noptloops
        self.kerneldatabase  = kerneldatabase       
        self.minhlambda      = minhlambda
        self.maxhlambda      = maxhlambda
        self.deltahlambda    = deltahlambda
        self.asconcentration = asconcentration

        if outputfilename is not None:
            self.outputfilename = outputfilename
        else:
            self.outputfilename = f"{model.name}.gpkde.out"

        self.parent.add_package(self)


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
        f = open(self.fn_path, "w")

        # Output filename
        f.write(f"{self.outputfilename}\n")

        # Domain origin
        for idb, b in enumerate(self.domainorigin):
            if idb == len(self.domainorigin)-1: 
                f.write(f"{b:10f}\n")
            else:
                f.write(f"{b:10f}\t")

        # Domain sizes
        for idb, b in enumerate(self.domainsize):
            if idb == len(self.domainsize)-1: 
                f.write(f"{b:10f}\n")
            else:
                f.write(f"{b:10f}\t")

        # Bin sizes
        for idb, b in enumerate(self.binsize):
            if idb == len(self.binsize)-1: 
                f.write(f"{b:10f}\n")
            else:
                f.write(f"{b:10f}\t")

        # Number of optimization loops 
        f.write(f"{self.noptloops:10d}\n")

        # Reconstruction with kernel database or brute force
        # 1: kerneldatabase
        # 0: bruteforcedensity
        if self.kerneldatabase:
            # Kernel database reconstruction
            f.write(f"1\n") # 1 for id into fortran

            # Database params: minh/lambda, daltah/lambda, maxh/lambda
            f.write(f"{self.minhlambda:10f}\t")
            f.write(f"{self.deltahlambda:10f}\t")
            f.write(f"{self.maxhlambda:10f}\n")

        else:
            # Brute force reconstruction
            f.write(f"0\n") # 0 for id into fortran
            
            # kernel params: minh/lambda, maxh/lambda
            f.write(f"{self.minhlambda:10f}\t")
            f.write(f"{self.maxhlambda:10f}\n")


        if self.asconcentration:
            f.write(f"1\n") # 1 for id into fortran
        else:
            f.write(f"0\n") # 0 for id into fortran


        # And close
        f.close()

        return
