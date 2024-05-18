'''
Configuration of MODPATH-RW sources
'''

# python
from collections import Counter
import numpy as np

# flopy
from flopy.pakbase import Package
from flopy.utils import Util3d, OptionBlock

# local 
from .utils import multipackage


class ModpathRWSrc( Package ):
    """
    MODPATH-RW Source Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class: ModpathRW) to which
        this package will be added.
    inputformat : str 
        The format of the specification for this source. Allowed values are:
          * AUX or AUXILIARY: build the source of particles from aux
            variables stored in the flow model. Configured ONLY via the 
            sources keyword argument (see below) and the ifaceoption.
          * SPEC or SPECIFIED: the user specifies time intervals,
            concentrations and other properties of the injection.
            Can be configured via the sources keyword or by 
            specifying the individual keyword arguments of the __init__ 
            function. 
    sources : list, tuple
        The specification of sources to be used. This parameter adopts
        different forms depending on the inputformat.

        * AUX/AUXILIARY format expects the data structure:

            sources = [
                [ budgetname, auxvar, particlesmass, (nx,ny,nz), speciesid ],
                ...
            ]

          or:

            sources = [
                [
                    budgetname,
                    [
                        [ auxvar_1, particlesmass_1, (nx,ny,nz)_1, speciesid_1 ],
                        [ auxvar_2, particlesmass_2, (nx,ny,nz)_2, speciesid_2 ],
                        ...
                    ]
                ],
                ...
            ]

          where:

            * budgetname    : str : The budget header/package from where to read the auxvar.
                                    The auxvar is interpreted as concentration, and the program 
                                    will extract the flow-rates from the budget. These two, and 
                                    the particlesmass, are combined to create a solute mass flux
                                    (a flux of particles). For mf6, this value can be the specific
                                    name given to the package stored in TXT2ID2 (see mf6io.docs). 
            * auxvar        : str : is the name of auxiliary variable.
            * particlesmass : float : the mass assigned to particles. It should be 
                                      greater than zero.
            * (nx,ny,nz) : int : a template for particles release. It could be a 
                                 single int for a uniform template. It should be 
                                 greater than one.
            * speciesid  : int : zero-based index of species to which the auxiliary 
                                 variable will be related. It will only be written to 
                                 the file if particlesmassoption == 2. It should be 
                                 greater or equal to zero.
          notes: 
            * The class will verify that pairs (budgetname,auxvar) are unique.
            * The minimum required specification is the budgetname and auxvar,
              and the rest of parameters are interpreted incrementally, meaning for
              example that if particlesmass is not given, none of the following 
              parameters will be interpreted and default values are assigned.

        * While using the SPEC/SPECIFIED inputformat, sources could be given as a list 
          of dictionaries with keys following the keyword arguments of the __init__ 
          function, starting from budgetname until particlesmass, like:
            
            sources = [
                {
                    'budgetname': 'WEL', 
                    'timeintervals' : [
                        [ts1,te1],
                        [ts2,te2],
                        ...
                    ],
                    'cellinput': 0, 
                    ...
                },
                ...
            ]

          notes: 
            * A detailed description of each keyword is given below.
            * if the format is SPEC/SPECIFIED and sources are given as a list of 
              dictionaries, then keyword arguments given to the constructor are 
              overriden. If sources=None, then a dictionary will be filled 
              with values obtained from the keyword arguments, creating a single 
              source specification.

    budgetname : str
        The budget header/package from where to extract flow-rates. For mf6, this value
        can be the specific name given to the package, stored in TXT2ID2 (see mf6io.docs). 
    timeintervals : list[ [float,float] ]
        A list of time intervals delimited by [tstart,tend]. Each interval defines 
        the duration of an injection with a given concentration. 
    cellinput : int
        Determines how to define the cells related to the source. Allowed values are:
          * 0: extract the cells from the budget file. Use all the cells related
               to the budgetname.
          * 1: cells are given as a list of cell ids in the cells keyword argument.
          * 2: cells are given as an array with u3d where 1 indicate that the cells
               should be considered.
        notes:
          * all of the given cells, regardless the specification format, will 
            release particles during the time intervals in which the flow-rate
            stored in budget is positive (into the aquifer).
    cells : list, np.ndarray
        The list of cell ids to be related to the injection. If cellinput == 2 then
        this should be the modelgrid array with 1 for cells to be considered 
        and 0 for those to be excluded.
    structured : bool
        The format to interpret the cells if cellinput == 1. True stands 
        for (lay row col) and False for the linear cellnumber.
    ifaceoption : int
        Flag to indicate whether an iface shall be applied for injection cells.
        iface can be specified individually for each cell.
          * 0: do not read iface
          * 1: read iface
        notes: 
          * if the source format is AUX/AUXILIARY, this flag will indicate 
            whether to extract or not the value of iface from the budget file.
            In the context of the source package, the iface value modifies
            a given particle template to be consistent with the idea of 
            releasing from a given cell face, compressing the corresponding
            orthogonal dimension.
    defaultiface : int
        Default value of iface for cells without an explicitly given value, 
        when ifaceoption==1. Only applies for the format SPEC/SPECIFIED. 
        A value of 0 preserves the given template as it is. 
    concpercell : int
        Flag to indicate that different concentrations are given for each 
        cell in the source. Only applies if cellinput == 1.
          * 0: all cells with the same concentration.
          * 1: cells with specific concentration.
    nspecies : int
        The number of species being injected. Defines how many concentrations 
        per time interval shall be written to the package file. At least 1.
    speciesid : list[int]
        The zero based indexes of species to which relate concentrations. Only 
        written to the package file if particlesmassoption == 2 (see ModpathRWSim).
    concentration : list
        The concentrations for each time interval. The list should have as 
        many entries as time intervals, and each entry is expected to as many 
        values as nspecies. For the specific case that concpercell == 1, 
        then each time interval should contain nspecies*ncells concentrations,
        with the rationale of nspecies index moving first.
    template : list[tuple(int)]
        The distribution of particles on a cell for a given species. 
        If templateoption==1, then one template is expected for different 
        species. These, in combination with the particlesmass, determine 
        the number of particles and overall the resolution 
        with which the source injection is characterized. 
    templateoption : int
        Flag to indicate whether one template is employed for all 
        species (templateoption = 0) or different templates for 
        each species (templateoption=1). 
    particlesmass : list[float]
        The particlemass for each species.
    nparticles : int
        If a template is not given, it will be filled with (nparticles nparticles nparticles).
    drape : int 
        Determines the treatment of particles placed on initially dry cells. 
        Allowed values are: 
          * 0 : particles are unreleased. 
          * 1 : particles are vertically transfered to the next uppermost active cell. 
        notes:
          The drape option only applies for MODFLOW models solved with the standard formulation. 
          Models solved with Newton-Raphson are not affected by this parameter.
    stringid : str, optional
        An id for the dispersion specification. If not given is automatically 
        filled with the package ftype and an integer id. 
    extension : str, optional
        File extension (default is 'src').
    """



    # Default options
    defaultcellinput      = 0
    defaultstructured     = True
    defaultifaceoption    = 0
    defaultdefaultiface   = 0
    defaultconcpercell    = 0
    defaultnspecies       = 1
    defaultspeciesid      = 0
    defaulttemplate       = [(2,2,2)]
    defaulttemplateoption = 0
    defaultparticlesmass  = 1.0
    defaultnparticles     = 2
    defaultdrape          = 0


    @staticmethod
    def _ftype():
        return 'SRC'


    @multipackage
    def __init__(
        self,
        model,
        inputformat    = 'aux',
        sources        = None,
        budgetname     = None, 
        timeintervals  = None,
        cellinput      = defaultcellinput,
        cells          = None, 
        structured     = defaultstructured, 
        ifaceoption    = defaultifaceoption, 
        defaultiface   = defaultdefaultiface, 
        concpercell    = defaultconcpercell, 
        nspecies       = defaultnspecies, 
        speciesid      = defaultspeciesid, 
        concentration  = None,
        template       = defaulttemplate, 
        templateoption = defaulttemplateoption, 
        particlesmass  = defaultparticlesmass, 
        nparticles     = defaultnparticles, 
        drape          = defaultdrape,
        stringid       = None, 
        extension      = 'src',
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

        # Save model data 
        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])
        self.shape3d = shape3d
        self.model        = model
        self.modelversion = self._parent.flowmodel.version 

        # Verify format type
        if ( not isinstance(inputformat,str) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type {str(type(inputformat))}"
                f" for the format specification. It should be str."
            )
        # Validate given format 
        if ( ( inputformat.upper() not in ['AUX','AUXILIARY','SPEC','SPECIFIED'] ) ):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid input format {str(inputformat)}."
                f" Allowed values are AUX (AUXILIARY) or SPEC (SPECIFIED)."
            )
        self.inputformat = inputformat.upper()

        # Save ifaceoption
        if ( ifaceoption is None ): 
            ifaceoption = 0
        if not isinstance( ifaceoption, int ): 
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for ifaceoption. It should be integer,"
                f" but {str(type(ifaceoption))} was given."
            )
        if ( ifaceoption not in [0,1] ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid ifaceoption {str(ifaceoption)}."
                f" The allowed values are 0 (do not read in budgets)"
                f" or 1 (read from budgets)."
            )
        self.ifaceoption = ifaceoption

        # Save drape
        if ( drape is None ): 
            drape = 0
        if not isinstance( drape, int ): 
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for drape. It should be integer,"
                f" but {str(type(drape))} was given."
            )
        if ( drape not in [0,1] ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid drape {str(drape)}."
                f" The allowed values are 0 or 1."
            )
        self.drape = drape


        # Process sources on the given format 
        # Will define self.sources
        if ( self.inputformat in ['AUX', 'AUXILIARY'] ):
            self._process_aux_format(sources=sources)
        elif ( self.inputformat in ['SPEC', 'SPECIFIED'] ):
            self._process_spec_format(
                    sources        = sources       ,
                    budgetname     = budgetname    ,
                    timeintervals  = timeintervals ,
                    cellinput      = cellinput     ,
                    cells          = cells         ,
                    structured     = structured    ,
                    ifaceoption    = ifaceoption   ,
                    defaultiface   = defaultiface  ,
                    concpercell    = concpercell   ,
                    nspecies       = nspecies      ,
                    speciesid      = speciesid     ,
                    concentration  = concentration ,
                    template       = template      ,
                    templateoption = templateoption,
                    particlesmass  = particlesmass ,
                    nparticles     = nparticles    ,
                    drape          = drape
                )

        # String id 
        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = ftype+str(self.id)

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

                # Write the source id
                f.write(f"{ins.stringid}\n")

                # Write format
                f.write(f"{ins.inputformat.upper()}\n")

                # Specifics for the aux input format
                if (
                    ( ins.inputformat.upper() == 'AUX' ) or 
                    ( ins.inputformat.upper() == 'AUXILIARY' )
                ):

                    # Create data format for AUX input 
                    fmts = []
                    fmts.append("{:20s}")  # auxvarname
                    fmts.append("{:.16f}") # particles mass
                    for nt in range(3):    # (nx,ny,nz)
                        fmts.append("{:3d}")
                    istart = 0
                    iend   = 5
                    if ( not hasattr(self.parent,'particlesmassoption') ):
                        raise Exception( 
                            f"{self.__class__.__name__}:"
                            f" The particlesmassoption was not found in parent package."
                            f" Did you define a ModpathRWSim package ?"
                        )
                    # Append solute id if particlesmassoption == 2
                    if ( 
                        ( self._parent.particlesmassoption == 2 ) 
                    ): 
                        fmts.append("{:4d}")
                        iend = 6
                    fmt  = " ".join(fmts) + "\n"
                    srcfmts = []
                    srcfmts.append("{:20s}") # srcname
                    srcfmt = " ".join(srcfmts) + "\n"


                    # Inform about the number of source budgets
                    f.write(f"{len(ins.uniquesources)}\n")

                    # (...) and loop over them
                    for isrc, src in enumerate( ins.uniquesources ):

                        # Write src name
                        if ( ins.ifaceoption == 1 ): 
                            if ( ins.drape == 0 ): 
                                f.write(f"{src.upper()}  {ins.ifaceoption}\n")
                            else:
                                f.write(f"{src.upper()}  {ins.ifaceoption}  {ins.drape}\n")
                        else:
                            if ( ins.drape == 0 ): 
                                f.write(f"{src.upper()}\n")
                            else:
                                f.write(f"{src.upper()}  0  {ins.drape}\n")


                        # Write the number of aux vars
                        f.write(f"{len(ins.allauxnames[isrc])}\n")
                        
                        # Write aux vars
                        for iaux, auxdata in enumerate(ins.alldataperaux[isrc]):
                            f.write(fmt.format(*auxdata[istart:iend]))


                # Specifics for the spec input format
                elif (
                    ( ins.inputformat.upper() == 'SPEC' ) or 
                    ( ins.inputformat.upper() == 'SPECIFIED' )
                ):
                    # Format for release template
                    tfmt = []
                    for nt in range(3):    # (nx,ny,nz)
                        tfmt.append("{:3d}")
                    tfmt = " " + " ".join(tfmt) + "\n"

                    # Inform about the number of source budgets
                    f.write(f"{len(ins.sources)}\n")

                    for isrc, src in enumerate( ins.sources ):
                    
                        if ( src['ifaceoption'] == 0 ):
                            if ( src['drape'] == 0 ):
                                # Only the budgetname
                                f.write(f"{src['budgetname'].upper()}\n")
                            else:
                                # Write with defaults 
                                f.write(f"{src['budgetname'].upper()}  0  0  {src['drape']}\n")
                        elif ( src['ifaceoption'] == 1 ):
                            if ( src['drape'] == 0 ):
                                # The budgetname plus iface options
                                f.write(
                                    f"{src['budgetname'].upper()}  1  {src['defaultiface']}\n"
                                )
                            else:
                                # Write all 
                                f.write(
                                    f"{src['budgetname'].upper()}  1  {src['defaultiface']}  {src['drape']}\n"
                                )

                        if ( src['cellinput'] in [0,2] ):
                            f.write(f"{src['cellinput']}\n")
                            # If zero, skip, will be read from budget internally by the program

                            # If 2, write with u3d
                            if (src['cellinput'] == 2 ):
                                src['cells'].get_file_entry()

                        elif ( src['cellinput'] == 1 ):

                            # If 1, write the list of cells
                            f.write(f"{src['cellinput']} {len(src['cells'])} {int(not src['structured'])} {src['concpercell']}\n")

                            if ( src['ifaceoption'] == 0 ):
                                # And the list
                                fmts = []
                                if src['structured']:
                                    fmts.append("{:9d}") # lay
                                    fmts.append("{:9d}") # row
                                    fmts.append("{:9d}") # col
                                else:
                                    fmts.append("{:9d}") # cellid
                                fmt = " " + " ".join(fmts) + "\n"
                                for oc in src['cells']:
                                    woc = np.array(oc).astype(np.int32)+1 # Correct the zero-based indexes
                                    f.write(fmt.format(*woc))

                            # Think what to do with specific ifaces when ifaceoption == 1
                            elif ( src['ifaceoption'] == 1 ):
                                # And the list
                                fmts = []
                                if src['structured']:
                                    fmts.append("{:9d}") # lay
                                    fmts.append("{:9d}") # row
                                    fmts.append("{:9d}") # col
                                else:
                                    fmts.append("{:9d}") # cellid
                                fmt = " " + " ".join(fmts) + "\n"
                                for oc in src['cells']:
                                    woc = np.array(oc).astype(np.int32)+1 # Correct the zero-based indexes
                                    f.write(fmt.format(*woc))

                        # Write nspecies and template option
                        f.write(f"{src['nspecies']}  {src['templateoption']}\n")

                        # Templates
                        # A single template for all species
                        if ( src['templateoption'] == 0 ):
                            f.write( tfmt.format(*src['template'][0]) )
                        # Species specific template
                        elif ( src['templateoption'] == 1 ):
                            for isp in range(src['nspecies']):
                                f.write( tfmt.format(*src['template'][isp]) )

                        # Format for particles mass
                        mfmt = []
                        for nm in range(src['nspecies']):   
                            mfmt.append("{:.16f}")
                        mfmt = " " + " ".join(mfmt) + "\n"

                        # Write the mass
                        f.write(mfmt.format(*src['particlesmass']))

                        # Give the species id if requested.
                        if ( not hasattr(self.parent,'particlesmassoption') ):
                            raise Exception( 
                                f"{self.__class__.__name__}:"
                                f" The particlesmassoption was not found in parent package."
                                f" Did you define a ModpathRWSim package ?"
                            )
                        if ( 
                            ( self._parent.particlesmassoption == 2 ) 
                        ): 
                            # Verify len
                            if ( src['nspecies'] != len(src['speciesid']) ): 
                                raise Exception(
                                    f"{self.__class__.__name__}:"
                                    f" The number of species is not consistent"
                                    f" with the number of given species ids."
                                    f" {str(src['nspecies'])} species were indicated"
                                    f" at nspecies and {str(len(src['speciesid']))}"
                                    f" were given as species ids." 
                            )

                            # The format
                            sfmt = []
                            for nm in range(src['nspecies']):   
                                sfmt.append("{:3d}")
                            sfmt = " " + " ".join(sfmt) + "\n"
                             
                            # and write 
                            f.write(sfmt.format(*src['speciesid']))

                        # And finally the time and concentration data
                        dfmt = []
                        dfmt.append("{:.16f}") # ts
                        dfmt.append("{:.16f}") # te
                        if ( src['cellinput'] in [0,2] ):
                            # One concentration column per species
                            for ns in range(src['nspecies']):   
                                dfmt.append("{:.16f}")
                        elif ( (src['cellinput'] == 1)and(src['concpercell']==0) ):
                            # One concentration column per species
                            for ns in range(src['nspecies']):   
                                dfmt.append("{:.16f}")
                        elif ( (src['cellinput'] == 1)and(src['concpercell']==1) ):
                            # nspecies concentrations per cell
                            for ns in range(src['nspecies']):   
                                dfmt.append("{:.16f}")
                        dfmt = " " + " ".join(dfmt) + "\n"

                        # Write the number of time intervals
                        f.write(f"{len(src['timeintervals'])}\n")

                        # And write the data
                        for itin, tin in enumerate(src['timeintervals']):
                            data = [tin[0], tin[1]]
                            data.extend(src['concentration'][itin])
                            f.write(dfmt.format(*data))

        # Done
        return


    def _process_aux_format( self, sources=None ):
        '''
        Specific protocols for interpretation of input format based on 
        auxiliary variables. 

        The function interprets the sources variable. 
        For the aux input format can have the structure

        sources = [
            [
                package_name, aux_var_name, particles_mass, (nx,ny,nz)
            ]
        ]

        or

        sources = [
            [
                package_name, 
                [
                    [ aux_var_name_1, particles_mass, (nx,ny,nz) ], 
                    [ aux_var_name_2, particles_mass, (nx,ny,nz) ], 
                ]
            ]
        ]

        Optionally, if the simulation is configured with particlesmassoption == 2, then 
        a species id can be given as a last parameter for the specification

        sources = [
            [
                package_name, aux_var_name, particles_mass, (nx,ny,nz), speciesid
            ]
        ]


        '''

        # Default variables for initialization of sources
        DEFAULTNP    = self.__class__.defaultnparticles
        DEFAULTMASS  = self.__class__.defaultparticlesmass
        DEFAULTSOLID = self.__class__.defaultspeciesid

        # Validate that sources were given
        if (
            ( sources is None )
        ):
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Input format AUX (AUXILIARY) requires the"
                f" specification of sources, but None was given."
            )

        # It should be a list or tuple
        if ( not isinstance( sources, (list, tuple) ) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Sources specification for input format"
                f" AUX (AUXILIARY) should be a list or tuple."
                f" {str(type(sources))} was given."
            )

        # Infer the number of sources
        NSOURCES = 0

        # If len sources is > 1, and the first entry is str, assume single source specification
        # If len sources == 1, verify that there is list or tuple inside and that the first entry is str
        lensources = len(sources)
        if ( lensources > 1 ):
            if isinstance( sources[0], str ):
                NSOURCES = 1
            else:
                # Iterate and verify that all entries are either list or tuple
                for s in sources:
                    if not isinstance( s, (list,tuple) ): 
                        raise ValueError(
                            f"{self.__class__.__name__}:"
                            f" Invalid sources specification for"
                            f" input format AUX (AUXILIARY). Verify structure."
                        )
                NSOURCES = lensources
        elif ( lensources == 1 ):
            if isinstance( sources[0], (list,tuple) ):
                if not isinstance( sources[0][0], str ): 
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Invalid sources specification for"
                        f" input format AUX (AUXILIARY). Verify structure."
                    )
                NSOURCES = 1
                sources  = sources[0]
            else:
                NSOURCES = 1

        # Short circuit just in case
        if ( NSOURCES == 0 ): 
            raise ValueError(
                f"{self.__class__.__name__}:"
                f" Invalid sources specification for input"
                f" format AUX (AUXILIARY). Empty sources list."
            )

        # Validate sources, verify auxiliary variables 
        sourceslist = []
        alldatalist = []

        # Read sources #
        for isrc in range(NSOURCES):

            # It NSOURCES is one, then source data is the whole spec
            if ( NSOURCES == 1): 
                sd = sources
            else:
                sd = sources[isrc]

            # Extract package name
            pname = sd[0]
            if not isinstance( pname, str ): 
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid type for package name in source"
                    f" specification. It should be str and"
                    f" {str(type(pname))} was given."
                )

            # With the package name, look for the package in the model
            pkg = self._parent.flowmodel.get_package(pname.lower())
            if pkg is None:
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Package name {pname.upper()} was not found in model."
                    f" Available packages are: {','.join(self._parent.flowmodel.get_package_list())}."
                )

            # Extract the list of aux vars for package
            # Note: auxlist filled only with str.upper()
            if self.modelversion == 'mf6':
                if not pkg.auxiliary.has_data():
                    raise Exception(
                        f"{self.__class__.__name__}:"
                        f" Auxiliary property without data for package {pkg.name[0]}."
                    )
                aux     = pkg.auxiliary.array[0]
                auxlist = [a.upper() for a in aux ]
                # Remove the aux identifiers
                auxlist = [ auxname for auxname in auxlist if auxname not in ('AUX','AUXILIARY')]
            else:
                # Extract possible aux names and preliminary checks for != MF6
                opts = pkg.options
                if isinstance(opts, OptionBlock):
                    auxlistcand = opts.auxillary
                elif isinstance(opts,(list,tuple)):
                    auxlistcand = pkg.options
                else:
                    raise TypeError(
                        f"{self.__class__.__name__}:"
                        f" Unexpected options type in package {pkg.name[0]}"
                        f" for the modelversion {self.modelversion}."
                    )
                # Verify at least one candidate
                if ( len(auxlistcand) == 0):
                    raise Exception(
                        f"{self.__class__.__name__}:"
                        f" No auxiliary variables were found in package {pkg.name[0]}."
                    )
                # Loop over aux candidates and create auxlist
                auxlist = []
                for ia, auxc in enumerate(auxlistcand):
                    # Clean the name ( str is assumed )

                    # Clean borders
                    cand = auxc.strip().upper()
                    
                    # If is the identifier str, go to the next
                    if ( cand in ['AUX','AUXILIARY'] ):
                        continue

                    # Break it
                    cand = cand.split()

                    # It it had a whitespace
                    if ( len(cand) == 2 ):
                        # Verify that the first word is the aux identifier, the second is the auxvar name
                        if ( cand[0].upper() in ['AUX','AUXILIARY'] ):
                            # Save it 
                            auxlist.append( cand[1].upper() )
                    elif ( len(cand) == 1):
                        # Save it 
                        auxlist.append( cand[0].upper() )
                    else:
                        # Probably the case in which auxvars have whitespace in between. 
                        raise Exception(
                            f"{self.__class__.__name__}:"
                            f" Unexpected definition of aux variable in {pkg.name[0]}."
                            f" Found name was {' '.join(cand)}."
                            f" Maybe multiple white spaces ?"
                        )

            # Short circuit just in case 
            if (len(auxlist) == 0):
                raise Exception(
                    f"{self.__class__.__name__}:"
                    f" No auxiliary variables were found,"
                    f" something wrong. Package {pkg.name[0]}."
                )
            
            # Now process the source specification #

            # Determine the format in which the aux variables were given 
            auxspec = sd[1]

            # If a single aux var spec
            if isinstance( auxspec, str ):
                # Initialize source variables with default values
                # List to store aux specs
                auxnames     = [auxspec.upper()]
                auxmasses    = [DEFAULTMASS]
                auxtemplates = DEFAULTNP*np.ones(shape=(1,3),dtype=np.int32)
                auxsolutes   = [DEFAULTSOLID]
                lensd        = len(sd)

                if ( (lensd >= 3 ) and (lensd <= 5 ) ):
                    # Verify particles mass
                    if not isinstance( sd[2], (int,float) ):
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid particle mass. It was expecting int/float, but"
                            f" {str(type(sd[2]))} was given."
                        )
                    if ( sd[2] <= 0 ): 
                        raise ValueError(
                            f"{self.__class__.__name__}:"
                            f" Invalid particle mass. It should be greater than zero." 
                            f" {str(sd[2])} was given."
                        )
                    # Save mass 
                    auxmasses[0] = sd[2]
                           
                    # Source specification as [auxvar,mass,template]
                    if lensd >= 4:
                        if isinstance(sd[3],(list,tuple)):
                            # Save template, if valid
                            for inp, npa in enumerate(sd[3]):
                                if ( npa <= 0 ): 
                                    raise ValueError(
                                        f"{self.__class__.__name__}:"
                                        f" Invalid particle template."
                                        f" All entries should be greater than zero."
                                        f" {str(npa)} was given at position {str(inp)}."
                                    )
                                auxtemplates[0,inp] = npa
                        elif isinstance( sd[3], int ):
                            # Assume uniform template
                            if ( sd[3] <= 0 ): 
                                raise ValueError(
                                    f"{self.__class__.__name__}:"
                                    f" Invalid particle template."
                                    f" All entries should be greater than zero."
                                    f" {str(sd[3])} was given."
                                )
                            auxtemplates[0,:] = sd[3]
                        else:
                            raise TypeError(
                                f"{self.__class__.__name__}:"
                                f" Invalid particles template."
                                f" It was expecting int/list/tuple, but"
                                f" {str(type(sd[3]))} was given."
                            )
                        # Source specification as [auxvar,mass,template,soluteid]
                        if lensd == 5:
                            if not isinstance(sd[4],int):
                                raise TypeError(
                                    f"{self.__class__.__name__}:"
                                    f" Invalid species id specification."
                                    f" It was expecting an integer, but"
                                    f" {str(type(sd[4]))} was given."
                                )
                            if ( sd[4] < 0 ):
                                raise ValueError(
                                    f"{self.__class__.__name__}:"
                                    f" Invalid species id specification."
                                    f" It was expecting a positive integer, but"
                                    f" {str(sd[4])} was given."
                                )
                            auxsolutes[0] = sd[4]
                else: 
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Invalid specification, exceeds maximum parameters length."
                    )

            # If a list of aux var specs
            elif isinstance( auxspec, (list,tuple)):
                # Initialize source variables with default values
                # List to store aux specs
                auxnames     = []
                auxmasses    = DEFAULTMASS*np.ones( shape=(len(auxspec),) )
                auxtemplates = DEFAULTNP*np.ones(shape=(len(auxspec),3),dtype=np.int32)
                auxsolutes   = DEFAULTSOLID*np.ones( shape=(len(auxspec),),dtype=np.int32 )

                # Loop over specifications and interpret the given values
                for iaspc, aspc in enumerate(auxspec):

                    if isinstance( aspc, (list,tuple) ):
                        auxspeclen = len(aspc)
                        if auxspeclen == 1:
                            # Source specification as [auxvar]
                            if not isinstance( aspc, str ):
                                raise TypeError(
                                    f"{self.__class__.__name__}:"
                                    f" Invalid auxiliary variable specification."
                                    f" It was expecting str, but"
                                    f" {str(type(aspc))} was given."
                                )
                            # Good, save it to the list of names
                            auxnames.append( aspc.upper() )
                        elif ( (auxspeclen >= 2 ) and (auxspeclen <= 4 ) ):
                            # Source specification as [auxvar,mass]
                            # Verify str for auxliary variable name 
                            if not isinstance( aspc[0], str ): 
                                raise TypeError(
                                    f"{self.__class__.__name__}:"
                                    f" Invalid auxiliary variable specification."
                                    f" It was expecting str, but"
                                    f" {str(type(aspc[0]))} was given."
                                )
                            # Good, save it to the list of names
                            auxnames.append( aspc[0].upper() )

                            # Verify particles mass
                            if not isinstance( aspc[1], (int,float) ):
                                raise TypeError(
                                    f"{self.__class__.__name__}:"
                                    f" Invalid particle mass. It was expecting"
                                    f" int/float, but"
                                    f" {str(type(aspc[1]))} was given."
                                )
                            if ( aspc[1] <= 0 ): 
                                raise ValueError(
                                    f"{self.__class__.__name__}:"
                                    f" Invalid particle mass. It should be greater than zero."
                                    f" {str(aspc[1])} was given."
                                )
                            # Good
                            auxmasses[iaspc] = aspc[1]
                           
                            # Continue reading 
                            if auxspeclen >= 3:
                                # Source specification as [auxvar,mass,template]
                                if isinstance(aspc[2],(list,tuple)):
                                    # Save template
                                    for inp, npa in enumerate(aspc[2]):
                                        if ( npa <= 0 ): 
                                            raise ValueError(
                                                f"{self.__class__.__name__}:"
                                                f" Invalid particle template."
                                                f" All entries should be greater than zero."
                                                f" {str(npa)} was given at position {str(inp)}."
                                            )
                                        auxtemplates[iaspc,inp] = npa
                                elif isinstance( aspc[2], int ):
                                    # Assume uniform template 
                                    if ( aspc[2] <= 0 ): 
                                        raise ValueError(
                                            f"{self.__class__.__name__}:"
                                            f" Invalid particle template."
                                            f" All entries should be greater than zero."
                                            f" {str(aspc[2])} was given."
                                        )
                                    auxtemplates[iaspc,:] = aspc[2]
                                else:
                                    raise TypeError(
                                        f"{self.__class__.__name__}:"
                                        f" Invalid particles template."
                                        f" It was expecting int/list/tuple, but"
                                        f" {str(type(aspc[2]))} was given."
                                    )
                                # Source specification as [auxvar,mass,template,soluteid]
                                if auxspeclen == 4:
                                    if not isinstance(aspc[3],int):
                                        raise TypeError(
                                            f"{self.__class__.__name__}:"
                                            f" Invalid species id specification."
                                            f" It was expecting an integer, but"
                                            f" {str(type(aspc[3]))} was given."
                                        )
                                    if ( aspc[3] < 0 ):
                                        raise ValueError(
                                            f"{self.__class__.__name__}:"
                                            f" Invalid species id specification."
                                            f" It was expecting a positive integer, but"
                                            f" {str(aspc[3])} was given."
                                        )
                                    auxsolutes[iaspc] = aspc[3]
                        else:
                            raise ValueError(
                                f"{self.__class__.__name__}:"
                                f" Invalid specification, maximum length is 4, but"
                                f" {str(len(aspc))} was given."
                            )
                    elif isinstance( aspc, str ):
                        # Is a str, save it to the list
                        auxnames.append( aspc.upper() )
                    else:
                        raise TypeError(
                            f"{self.__class__.__name__}:"
                            f" Invalid auxiliary variable"
                            f" specification. It was expecting str."
                            f" {str(type(aspc))} was given."
                        )
            else:
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid auxiliary variable specification."
                    f" It was expecting type str or a list/tuple."
                    f" {str(type(auxspec))} was given."
                )


            # Validate given auxnames
            for ian, an in enumerate(auxnames):
                acount = auxlist.count(an.upper())
                # If the aux name is more than once 
                if acount > 1:
                    raise Exception(
                        f"{self.__class__.__name__}:"
                        f" Auxiliary name {an.upper()} was given"
                        f" more than once in package"
                        f" {pname.upper()}. Check definition."
                    )
                # If the given name is not in the list of aux vars of the package
                if (an.upper() not in auxlist):
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Auxiliary name {an.upper()} was not found"
                        f" in aux variables of package"
                        f" {pname.upper()}. Available variables are: {','.join(auxlist)}."
                    )
                # Definition is consistent, save the pair [sourcename,auxname]
                tempair  = [pname.upper(),an.upper()]
                sourceslist.append(tempair)
                tempdata = [
                    an.upper(),
                    auxmasses[ian],
                    auxtemplates[ian,0],
                    auxtemplates[ian,1],
                    auxtemplates[ian,2],
                    auxsolutes[ian] + 1 # Notice the + 1, fortran index !
                ]
                alldatalist.append(tempdata)
        # end read sources # 

        # Validate that combinations of [sourcename, auxname] are unique 
        counter       =  Counter(frozenset(sl) for sl in sourceslist )
        counterstatus = [counter[frozenset(sl)] == 1 for sl in sourceslist] # True or False

        # Break if some a combination was given more than once
        if not np.all( counterstatus ):
            raise Exception(
                f"{self.__class__.__name__}:"
                f" Some pairs [package_name,aux_var_name] are repeated."
                f" Combination should be unique for all source packages."
            )

        # Combination of sourcenames,auxnames is valid
        # Extract unique source names
        self.uniquesources = np.unique( np.array(sourceslist)[:,0] )
        self.allauxnames   = []
        self.alldataperaux = []
        for ius, us in enumerate( self.uniquesources ):
            auxnamespersource = []
            dataperaux        = []
            for isl, sl in enumerate(sourceslist):
                if sl[0] == us:
                    auxnamespersource.append(sl[1])
                    data = np.array(alldatalist[isl],dtype=object)
                    dataperaux.append(data)
            # Indexed as uniquesources
            self.allauxnames.append( auxnamespersource )
            self.alldataperaux.append( dataperaux )
        
        # If mf6, also get the source kind for 
        # the given srcnames, just in case
        self.sourcestype = []
        if self.modelversion == 'mf6':
            for us in self.uniquesources:
                pkg = self._parent.flowmodel.get_package(pname.lower())
                self.sourcestype.append( pkg.package_type.upper() )

        # Save sources spec
        self.sources = sources

        # And done
        return


    def _process_spec_format( self, *args, **kwargs ):
        '''
        Specific protocols for interpretation of input format based on 
        user specified variables

        Will define self.sources with a variable following the format:
                    
            sources = [
                    {
                        'budgetname'    : 'WEL-1',
                        'timeintervals' : [
                            [ 0, 1*365*86400 ],
                        ],
                        'cellinput'     : 1,
                        'cells'         : [48],
                        'structured'    : False,
                        'ifaceoption'   : 0, 
                        'defaultiface'  : 0,
                        'concpercell'   : 0,
                        'nspecies'      : 1,
                        'speciesid'     : [0],
                        'templateoption': 0,
                        'template'      : [(2,2,2)],
                        'particlesmass' : [1.0],
                        'concentration' : [
                                [55],
                        ],
                    },
                ]

        '''
      
        # If sources is None, define one with 
        # the keyword arguments 
        if ( kwargs['sources'] is None ):
            # Just to be consisten, remove the source kw
            kwargs.pop( 'sources', None )
            sources = [ kwargs ]
        else:
            sources = kwargs['sources']

        # Do some validation for the given sources
        # This format is more rigid, force reading specs from dicts
        if ( not isinstance(sources,(list,dict)) ):
            raise TypeError(
                f"{self.__class__.__name__}:"
                f" Invalid type for sources in specified input format,"
                f" expected list or dict. {str(type(sources))} was given."
            )
        if ( isinstance(sources,(dict)) ):
            sources = [sources]
        if ( isinstance(sources,list) ):
            for src in sources:
                if ( not isinstance( src, dict ) ): 
                    raise TypeError(
                        f"{self.__class__.__name__}:"
                        f" Invalid type for source in specified input format, expected dict."
                        f" {str(type(src))} was given."
                    )

        # Do more advanced validation for each source specification 
        # Consider recycling the machinery for validating params in mf6 classes
        for src in sources:

            # Get the current source keys
            keys = src.keys()

            # Verify keys specification

            # budgetname
            if ('budgetname' not in keys ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The key"
                    f" budgetname was not found in specification."
                )
            if ( src['budgetname'] is None ): 
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The key"
                    f" budgetname is None."
                )
            if ( not isinstance( src['budgetname'], str ) ): 
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The key"
                    f" budgetname is not a str."
                    f" {str(type(src['budgetname']))} was given."
                )
            # It could be considered a check versus the packages available at the model 
            # but headers could have different names e.g. "CONSTANT HEAD" in mf2005

            # nspecies
            if ('nspecies' not in keys ):
                # Give the defaults
                src['nspecies'] = self.__class__.defaultnspecies
                if ('speciesid' not in keys ):
                    src['speciesid'] = self.__class__.defaultspeciesid
            if ( not isinstance( src['nspecies'], int ) ): 
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The key"
                    f" nspecies should be int. {str(type(src['nspecies']))} was given."
                )
            if ( src['nspecies'] <= 0 ): 
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The nspecies should"
                    f" be a positive integer."
                    f" {str(src['nspecies'])} was given."
                )
            nsp = src['nspecies']
            
            # speciesid
            # This is finally validated while writing, 
            # as is the only place where particlesmassoption 
            # can be accessed.
            if ( 'speciesid' not in keys ):
                # initialize with the same number of ids as the 
                # one given in nspecies. Still it might need validation 
                # with respect to the total number of given concentrations. 
                src['speciesid'] = [ sid for sid in range(src['nspecies']) ]
            if ( isinstance( src['speciesid'], int ) ): 
                # Create a list
                src['speciesid'] = [src['speciesid']]
            if ( isinstance( src['speciesid'], list ) ): 
                # Notice the + 1, fortran index !
                src['speciesid'] = np.array( src['speciesid' ] ).astype(np.int32) + 1
                src['speciesid'] = src['speciesid'].tolist()

            # particlesmass
            if ('particlesmass' not in keys ):
                # Give the default
                src['particlesmass'] = self.__class__.defaultparticlesmass
            if ( not isinstance( src['particlesmass'], list ) ):
                if ( not isinstance( src['particlesmass'], (int,float) ) ): 
                    raise TypeError(
                        f"{self.__class__.__name__}:"
                        f" Invalid source specification. The given"
                        f" particles mass is not int/float."
                        f" {str(type(src['particlesmass']))} was given."
                    )
                if ( nsp == 1 ): 
                    src['particlesmass'] = [src['particlesmass']]
                else:
                    # Assign the default multiple times ( or throw exception ? ) 
                    src['particlesmass'] = list(np.ones(shape=(nsp,))*src['particlesmass'])
            else:
                if( not ( len( src['particlesmass'] ) == nsp ) ):
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Invalid source specification. The number of"
                        f" elements in particlesmass list"
                        f" is not consistent with nspecies."
                        f" Nspecies is {str(nsp)} and"
                        f" {str(len(src['particlesmass']))} masses were given."
                    )

            # cellinput
            if ('cellinput' not in keys ):
                # Give the default
                src['cellinput'] = self.__class__.defaultcellinput
            if ( src['cellinput'] not in [0,1,2] ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The cellinput"
                    f" parameter is not 0, 1 or 2."
                    f" {str(src['cellinput'])} was given."
                )

            # structured 
            if ('structured' not in keys ):
                # Give the default
                src['structured'] = self.__class__.defaultstructured

            # ifaceoption
            if ('ifaceoption' not in keys ):
                # Give the default
                src['ifaceoption'] = self.__class__.defaultifaceoption
            if ( src['ifaceoption'] not in [0,1] ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The ifaceoption"
                    f" parameter is not 0 or 1."
                    f" {str(src['ifaceoption'])} was given."
                )

            # defaultiface
            if ('defaultiface' not in keys ):
                # Give the default
                src['defaultiface'] = self.__class__.defaultdefaultiface
            if ( not isinstance( src['defaultiface'], int ) ):
                raise TypeError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The defaultiface is"
                    f" not an integer."
                    f" {str(type(src['defaultiface']))} was given."
                )
            if ( not (0<= src['defaultiface']<=6) ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. Allowed values"
                    f" for defaultiface are from 0 to 6."
                    f" {str(src['defaultiface'])} was given."
                )

            # concpercell
            if ('concpercell' not in keys ):
                # Give the default
                src['concpercell'] = self.__class__.defaultconcpercell
            if ( src['concpercell'] not in [0,1] ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The concpercell"
                    f" parameter is not 0 or 1."
                    f" {str(src['concpercell'])} was given."
                )
            concpercell = src['concpercell']

            # drape
            if ('drape' not in keys ):
                # Give the default
                src['drape'] = self.__class__.defaultdrape
            if ( src['drape'] is None ):
                src['drape'] = self.__class__.defaultdrape
            if ( src['drape'] not in [0,1] ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The drape"
                    f" parameter is not 0 or 1."
                    f" {str(src['drape'])} was given."
                )


            # templateoption 
            if ('templateoption' not in keys ):
                # Give the default
                src['templateoption'] = self.__class__.defaulttemplateoption
            if ( src['templateoption'] not in [0,1] ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The templateoption"
                    f" parameter is not 0 or 1."
                    f" {str(src['templateoption'])} was given."
                )

            # template
            if ('template' not in keys ):
                # Give the default
                src['template'] = self.__class__.defaulttemplate
            # Requires some default definition and validation

            # cells
            if ( ( src['cells'] is None ) and ( src['cellinput'] in [1,2] ) ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The cell input option"
                    f" {str(src['cellinput'])} requires cells to be"
                    f" specified but None was given."
                )
            if ( src['cellinput'] == 0 ):
                # Same concentration for all cells
                ncellsforinput = 1
            if ( src['cellinput'] == 1 ):
                if ( not isinstance( src['cells'], list ) ):
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Invalid source specification. The cell"
                        f" input option is 1 but no list of cells was given."
                    )
                # Requires something to validate versus the structured parameter
                cells          = np.array(src['cells'])
                ncellsforinput = cells.shape[0]
                # If specifying the same concentration for all cells 
                # then update ncellsforinput to only one cell 
                if ( src['concpercell'] == 0 ):
                    ncellsforinput = 1
                
            if ( src['cellinput'] == 2 ):
                # Same concentration for all cells
                cells = src['cells']
                src['cells'] = Util3d(
                    self.model,
                    self.shape3d,
                    np.int32,
                    cells,
                    name="SRCCELLS",
                    locat=self._parent.multipackage[ftype]['unitnumber'], 
                )
                ncellsforinput = 1

            # timeintervals
            if ('timeintervals' not in keys ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The key"
                    f" timeintervals was not found in specification."
                )
            if ( not isinstance( src['timeintervals'], list ) ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The key timeintervals should be"
                    f" a list, but {str(type(src['timeintervals']))} was given."
                )
            else:
                tintervals = np.array( src['timeintervals'] )
                # Handle the case when given only one tstart, tend
                if ( len( tintervals.shape ) == 1 ):
                    tintervals = np.expand_dims( tintervals, axis=0 )
                if ( not ( tintervals.shape[1] == 2 ) ): 
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Invalid source specification. The key"
                        f" timeintervals should be have a"
                        f" shape of two in axis=1 indicating tstart, tend."
                    )
                ntin = tintervals.shape[0]
        
                # Give it back
                src['timeintervals'] = tintervals.tolist()


            # concentration
            if ('concentration' not in keys ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The key"
                    f" concentration was not found in specification."
                )
            if ( not isinstance( src['concentration'], list ) ):
                raise ValueError(
                    f"{self.__class__.__name__}:"
                    f" Invalid source specification. The key concentration should be"
                    f" a list, but {str(type(src['concentration']))} was given."
                )
            else:
                # To numpy array 
                conc = np.array(src['concentration'])
                if ( len( conc.shape ) == 1 ): 
                    conc = np.expand_dims( conc, axis=1 )
                if ( conc.shape[0] != ntin ): 
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Invalid source specification. The number of given concentrations"
                        f" ({str(conc.shape[0])}) is not consistent with the"
                        f" number of time intervals ({str(ntin)})."
                    )
                # Verifies how many concentrations 
                if ( conc.shape[1] != nsp*ncellsforinput ):
                    raise ValueError(
                        f"{self.__class__.__name__}:"
                        f" Invalid source specification. The given"
                        f" number of concentrations/columns ({str(conc.shape[1])}) is not consistent"
                        f" with the expected number ({str(nsp*ncellsforinput)})."
                        f" The expected number of concentrations is obtained as"
                        f" the product between the given number of species ({str(nsp)})"
                        f" and the expected number of cell-specific concentrations ({str(ncellsforinput)}),"
                        f" due to the given concpercell ({str(concpercell)}) parameter."
                    )
                # Pass it back 
                src['concentration'] = conc.tolist()


        # Save to sources
        self.sources = sources

        # And done
        return
