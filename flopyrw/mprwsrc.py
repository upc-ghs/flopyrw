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
from .utils import count_instances # Increment COUNTER


class ModpathRWSrc( Package ):
    """
    MODPATH-RW Source Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class: ModpathRW) to which
        this package will be added.
    format : str 
        The format specification for this source. Allowed values are:
            * AUX or AUXILIARY: defines the source from aux variables stored in the flow model. 
                                Configured ONLY via the sources keyword argument and ifaceoption.
            * SPEC or SPECIFIED: user specifies time intervals, concentrations and other properties
                                 of the injection. Can be configured via the sources keyword or by 
                                 specifying the individual keyword arguments of the __init__ function. 
    sources : list, tuple
        The specification of sources to be used. 
        * While using the AUX/AUXILIARY format parameter expects the data structure:

                sources = [
                    [ package_name, aux_var_name, particles_mass, (nx,ny,nz), speciesid ],
                    ...
                ]

            or:

                sources = [
                    [
                        package_name,
                        [
                            [ aux_var_name_1, particles_mass_1, (nx,ny,nz)_1, speciesid_1 ],
                            [ aux_var_name_2, particles_mass_2, (nx,ny,nz)_2, speciesid_2 ],
                            ...
                        ]
                    ],
                    ...
                ]

            where:

                * package name : str : indicates the package name from where to extract the aux variables
                * aux_var_name : str : is the name of auxiliary variable
                * particles_mass : float : the mass assigned to particles. It should be greater than zero.
                * (nx,ny,nz) : int : a template for particles release. It could be a single int
                                     for a uniform template. It should be greater than one.
                * speciesid : int : zero-based index of species to which the auxiliary variable will be related.
                                    It will only be written to the file if particlesmassoption == 2.
                                    It should be greater or equal to zero.
            notes: 
                * The class will verify that pairs (package_name,aux_var_name) are unique.
                * The minimum required specification is the package_name and aux_var_name, and 
                  the rest of parameters are interpreted incrementally, meaning for example that
                  if particles_mass is not given, none of the following will be interpreted 
                  and default values are assigned
        * While using the SPEC/SPECIFIED input format, sources could be given as a list of dictionaries 
          with keys as the keyword parameters of the __init__ function starting from budgetname until
          particlesmass, like:
            
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
            notes: if the format is SPEC/SPECIFIED and sources is given as before, then the keyword 
                   parameters are overriden. If sources=None, then a dictionary will be filled by the 
                   class with the values obtained from the keyword arguments, creating a single source.

    budgetname : str
        The budget header from where to extract flow-rates. For mf6 is the name given to the package, 
        stored in TXTID2 (see mf6io.docs). 
    timeintervals : list[ [float,float] ]
        A list of time intervals delimited by [tstart,tend]. Each interval defines the duration of 
        an injection at a given concentration. 
    cellinput : int
        Determines how to define the cells related to the source. Allowed values are:
            * 0: extract the cells from the budget file. Use all the cells related to the budgetname.
            * 1: cells are given as a list of cell ids in the cells keyword argument.
            * 2: cells are given as an array with u3d where 1 indicate that the cells should be considered.
    cells : list, np.array
        The list of cell ids to be related to the injection. If cellinput==2 then this should be the modelgrid 
        array with 1 for cells to be considered and 0 for those to be excluded.
    structured : bool
        The format to interpret the cells if cellinput == 1. True stands for (lay row col) and False 
        for the linear cellnumber
    ifaceoption : int
        Flag to indicate whether an iface shall be applied for injection cells. iface can be specified 
        individually for each cell.
            * 0: do not read iface
            * 1: read iface
        notes: if the source format is AUX/AUXILIARY, this flag will indicate whether to extract or not
               the value of iface from the budget file. In the context of the source package, the iface
               value modifies a given particle template to be consistent with the idea of releasing 
               from a given cell face, compressing the corresponding orthogonal dimension.
    defaultiface : int
        Default value of iface for cells without an explicitly given value, when ifaceoption==1. Only 
        applies for the format SPEC/SPECIFIED. A value of 0 preserves the given template as it is. 
    concpercell : int
        Flag to indicate that different concentrations are given for each cell in the source. Only
        applies if cellinput == 1.
            * 0: all cells with the same concentration
            * 1: cells with specific concentration
    nspecies : int
        The number of species being injected. Defines how many concentrations per time interval shall 
        be written to the package file. At least 1.
    speciesid : list[int]
        The zero based indexes of species to which relate concentrations. Only written to the package
        file if particlesmassoption==2.
    concentration : list
        The concentrations for each time interval. The list should have as many entries as time intervals, 
        and each entry is expected to as many values as nspecies. For the specific case that concpercell == 1, 
        then each time interval should contain nspecies*ncells concentrations, with the rationale of nspecies
        index moving first.
    template : list[tuple(int)]
        The distribution of particles on a cell for a given species. If templateoption==1,
        then one template is expected for different species. These, in combination with 
        the particlesmass, determine the number of particles and overall the resolution 
        with which the source injection is characterized. 
    templateoption : int
        Flag to indicate whether one template is employed for all species (templateoption = 0)
        or different templates for each species (templateoption=1). 
    particlesmass : list[float]
        The particlemass for each species.
    nparticles : int
        If a template is not given, it will be filled with (nparticles nparticles nparticles)
    id : int 
        A positive integer identifier. If None given is automatically assigned with 
        the instance counter. 
    stringid : str, optional
        An id for the dispersion specification. If not given is automatically filled with 
        DSP+str(id), where id is an internal use counter. 
    extension : str, optional
        File extension (default is 'src').
    """


    # Class properties
    COUNTER    = 0
    UNITNUMBER = 0
    INSTANCES  = []


    @count_instances
    def __init__(
        self,
        model,
        format         = 'aux',
        sources        = None,
        budgetname     = None, 
        timeintervals  = None,
        cellinput      = 0,
        cells          = None, 
        structured     = True, 
        ifaceoption    = 0,
        defaultiface   = 0,
        concpercell    = 0,
        nspecies       = 1,
        speciesid      = 0,
        concentration  = None,
        template       = [(2,2,2)],
        templateoption = 0, 
        particlesmass  = 1.0,
        nparticles     = 2,
        id             = None, 
        stringid       = None, 
        extension      = 'src',
    ):
       
        # Define UNITNUMBER if the first instance is created
        if self.__class__.COUNTER == 1:
            unitnumber = model.next_unit()
            super().__init__(model, extension, "SRC", unitnumber)
            self.__class__.UNITNUMBER = self.unit_number[0]
        else:
            self._parent = self.INSTANCES[0]._parent

        # Save model data 
        self.model        = model
        self.modelversion = self._parent.flowmodel.version 

        # Verify format type
        if ( not isinstance(format,str) ):
            raise TypeError(
                self.__class__.__name__ + ':' + 
                ' Invalid type ' +str(type(format))+ ' for the format specification.' +
                '. It should be str. '
            )
        # Validate given format 
        if ( ( format.upper() not in ['AUX','AUXILIARY','SPEC','SPECIFIED'] ) ):
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid input format ' +str(format)+
                '. Allowed values are AUX (AUXILIARY) or SPEC (SPECIFIED). '
            )
        self.format = format.upper()
        # Save ifaceoption
        if ( ifaceoption not in [0,1] ): 
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Invalid ifaceoption ' +str(ifaceoption)+
                '. The allowed values are 0 (do not read in budgets) or 1 (read from budgets).'
            )
        self.ifaceoption = ifaceoption


        # Process sources on the given format 
        # Will define self.sources
        if ( self.format in ['AUX', 'AUXILIARY'] ):
            self._process_aux_format(sources=sources)
        elif ( self.format in ['SPEC', 'SPECIFIED'] ):
            self._process_spec_format(sources=sources)

        # Define id
        if (id is not None): 
            self.id = id
        else:
            self.id = self.__class__.COUNTER

        if (stringid is not None): 
            self.stringid = stringid
        else:
            self.stringid = 'SRC'+str(self.__class__.COUNTER)


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


        # Create data format for AUX input 
        fmts = []
        fmts.append("{:20s}")  # auxvarname
        fmts.append("{:.6f}")  # particles mass
        for nt in range(3):    # (nx,ny,nz)
            fmts.append("{:3d}")
        istart = 0
        iend   = 5
        # Append solute id is particlesmassoption == 2
        if ( 
            ( self.INSTANCES[0]._parent.particlesmassoption == 2 ) 
        ): 
            fmts.append("{:4d}")
            iend = 6
        fmt  = " ".join(fmts) + "\n"
        srcfmts = []
        srcfmts.append("{:20s}") # srcname
        srcfmt = " ".join(srcfmts) + "\n"


        # Write how many INSTANCES 
        # and loop over them
        f.write(f"{self.COUNTER}\n")

        for ins in self.INSTANCES:

            # Write the source id
            f.write(f"{ins.stringid}\n")

            # Write format
            f.write(f"{ins.format.upper()}\n")

            # Specifics for the aux input format
            if (
                ( ins.format.upper() == 'AUX' ) or 
                ( ins.format.upper() == 'AUXILIARY' )
            ):

                # Inform about the number of source budgets
                f.write(f"{len(ins.uniquesources)}\n")

                # (...) and loop over them
                for isrc, src in enumerate( ins.uniquesources ):

                    # Write src name
                    if ( ins.ifaceoption == 1 ): 
                        f.write(f"{src.upper()}    {ins.ifaceoption}\n")
                    else:
                        f.write(f"{src.upper()}\n")

                    # Write the number of aux vars
                    f.write(f"{len(ins.allauxnames[isrc])}\n")
                    
                    # Write aux vars
                    for iaux, auxdata in enumerate(ins.alldataperaux[isrc]):
                        f.write(fmt.format(*auxdata[istart:iend]))


            # Specifics for the spec input format
            elif (
                ( ins.format.upper() == 'SPEC' ) or 
                ( ins.format.upper() == 'SPECIFIED' )
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
                        # Only the budgetname
                        f.write(f"{src['budgetname'].upper()}\n")
                    elif ( src['ifaceoption'] == 1 ):
                        # The budgetname plus iface options
                        f.write(f"{src['budgetname'].upper()}    1    {src['defaultiface']}\n")

                    if ( src['cellinput'] in [0,2] ):
                        f.write(f"{src['cellinput']}\n")
                        # If zero, skip, will be read from budget internally by the program

                        # If 2, write with u3d
                        if (src['cellinput'] == 2 ):
                            src['cells'].get_file_entry()

                    elif ( src['cellinput'] == 1 ):

                        # If 1, write the list of cells
                        f.write(f"{src['cellinput']}  {len(src['cells'])}  {int(not src['structured'])}  {src['concpercell']}\n")

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
                                if src['structured']:
                                    f.write(fmt.format(*woc))
                                else:
                                    f.write(fmt.format(woc))

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
                                if src['structured']:
                                    f.write(fmt.format(*woc))
                                else:
                                    f.write(fmt.format(woc))

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
                        mfmt.append("{:.6f}")
                    mfmt = " " + " ".join(mfmt) + "\n"

                    # Write the mass
                    f.write(mfmt.format(*src['particlesmass']))

                    # Give the species id if requested. 
                    if ( 
                        ( self.INSTANCES[0]._parent.particlesmassoption == 2 ) 
                    ): 
                        # The format
                        sfmt = []
                        for nm in range(src['nspecies']):   
                            sfmt.append("{:3d}")
                        sfmt = " " + " ".join(sfmt) + "\n"
                        # and write 
                        f.write(sfmt.format(*src['speciesid']))

                    # And finally the time and concentration data
                    dfmt = []
                    dfmt.append("{:.6f}") # ts
                    dfmt.append("{:.6f}") # te
                    if ( src['cellinput'] in [0,2] ):
                        # One concentration column per species
                        for ns in range(src['nspecies']):   
                            dfmt.append("{:.6f}")
                    elif ( (src['cellinput'] == 1)and(src['concpercell']==0) ):
                        # One concentration column per species
                        for ns in range(src['nspecies']):   
                            dfmt.append("{:.6f}")
                    elif ( (src['cellinput'] == 1)and(src['concpercell']==1) ):
                        # nspecies concentrations per cell
                        for ns in range(src['nspecies']):   
                            dfmt.append("{:.6f}")
                    dfmt = " " + " ".join(dfmt) + "\n"

                    # Write the number of time intervals
                    f.write(f"{len(src['timeintervals'])}\n")

                    # And write the data
                    for itin, tin in enumerate(src['timeintervals']):
                        data = [tin[0], tin[1]]
                        data.extend(src['concentration'][itin])
                        f.write(dfmt.format(*data))

                    # And done

        # And close
        f.close()


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
        DEFAULTNP    = 2
        DEFAULTMASS  = 1.0
        DEFAULTSOLID = int(0)

        # Validate that sources were given
        if (
            ( sources is None )
        ):
            raise ValueError(
                self.__class__.__name__ + ':' + 
                ' Input format AUX (AUXILIARY) requires the specification of sources. '+
                'None was given.' 
            )

        # It should be a list or tuple
        if ( not isinstance( sources, (list, tuple) ) ):
            raise TypeError(
                self.__class__.__name__ + ':' + 
                ' Sources specification for input format AUX (AUXILIARY) should be a list or tuple. ' + 
                str(type(sources)) + ' was given.'
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
                            self.__class__.__name__ + ':' + 
                            ' Invalid sources specification for input format AUX (AUXILIARY). '+
                            'Verify structure.' 
                        )
                NSOURCES = lensources
        elif ( lensources == 1 ):
            if isinstance( sources[0], (list,tuple) ):
                if not isinstance( sources[0][0], str ): 
                    raise ValueError(
                        self.__class__.__name__ + ':' + 
                            ' Invalid sources specification for input format AUX (AUXILIARY). '+
                            'Verify structure.' 
                    )
                NSOURCES = 1
                sources  = sources[0]
            else:
                NSOURCES = 1

        # Short circuit just in case
        if ( NSOURCES == 0 ): 
            raise ValueError(
                self.__class__.__name__ + ':' + 
                    ' Invalid sources specification for input format AUX (AUXILIARY). '+
                    'Empty sources list.' 
            )

        # Validate sources, verify auxiliary variables 
        sourceslist = []
        alldatalist = []

        ## read sources ##
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
                    self.__class__.__name__ + ':' + 
                        ' Invalid type for package name in source specification. ' + 
                        'It should be str and ' +str(type(pname)) + ' was given.' 
                )

            # With the package name, look for the package in the model
            pkg = self._parent.flowmodel.get_package(pname.lower())
            if pkg is None:
                raise Exception(
                    self.__class__.__name__ + ':' + 
                    ' Package name ' + pname.upper() + ' was not found in model. ' 
                    'Available packages are: ' +  ','.join(self._parent.flowmodel.get_package_list()) 
                )

            # Extract the auxlist for package
            # Note: auxlist filled only with str.upper()
            # After this block the list auxlist is expected to be defined
            if self.modelversion == 'mf6':
                if not pkg.auxiliary.has_data():
                    raise Exception(
                        self.__class__.__name__ + ':' + 
                        ' Auxiliary property without data for package ' + pkg.name[0] 
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
                        self.__class__.__name__ + ':' + 
                        ' Unexpected options type in package ' + pkg.name[0] + 
                        'for the modelversion '+ self.modelversion
                    )
                # Veriy at least one candidate
                if ( len(auxlistcand) == 0):
                    raise Exception(
                        self.__class__.__name__ + ':' + 
                        ' No auxiliary variables were found in package ' + pkg.name[0] 
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
                            self.__class__.__name__ + ':' + 
                            ' Unexpected definition of aux variable in ' + pkg.name[0] +
                            'Found name was ' + ' '.join(cand) + '. Maybe multiple white spaces ?' 
                        )
            # Short circuit just in case 
            if (len(auxlist) == 0):
                raise Exception(
                    self.__class__.__name__ + ':' + 
                    ' No auxiliary variables were found, something wrong. Package ' + pkg.name[0]
                )
            
            # Now process the sources specification

            # Determine the format in which the aux variables were given 
            auxspec = sd[1]

            # If a single aux var spec
            if isinstance( auxspec, str ):
                # List to store aux specs
                auxnames     = [auxspec.upper()]
                # Initialize source variables with default values
                auxmasses    = [DEFAULTMASS]
                auxtemplates = DEFAULTNP*np.ones(shape=(1,3),dtype=np.int32)
                auxsolutes   = [DEFAULTSOLID]

                lensd = len(sd)
                if ( (lensd >= 3 ) and (lensd <= 5 ) ):
                    # Verify particles mass
                    if not isinstance( sd[2], (int,float) ):
                        raise TypeError(
                            self.__class__.__name__ + ':' + 
                            ' Invalid particle mass. It was expecting int/float, but ' + 
                            str(type(sd[2]))+ ' was given.'
                        )
                    if ( sd[2] <= 0 ): 
                        raise ValueError(
                            self.__class__.__name__ + ':' + 
                            ' Invalid particle mass. It should be greater than zero. ' +  
                            str(sd[2])+ ' was given.'
                        )
                    # Save mass 
                    auxmasses[0] = sd[2]
                           
                    # Continue reading 
                    if lensd >= 4:
                        # Source specification as [auxvar,mass,template]
                        if isinstance(sd[3],(list,tuple)):
                            # Save template
                            for inp, npa in enumerate(sd[3]):
                                auxtemplates[0,inp] = npa
                        elif isinstance( sd[3], int ):
                            # Assume uniform template 
                            auxtemplates[0,:] = sd[3]
                        else:
                            raise TypeError(
                                self.__class__.__name__ + ':' + 
                                ' Invalid particles template. It was expecting int/list/tuple, but ' + 
                                str(type(sd[3]))+ ' was given.'
                            )
                        # Source specification as [auxvar,mass,template,soluteid]
                        if lensd == 5:
                            if not isinstance(sd[4],int):
                                raise TypeError(
                                    self.__class__.__name__ + ':' + 
                                    ' Invalid species id specification. It was expecting an integer, but ' + 
                                    str(type(sd[4]))+ ' was given.'
                                )
                            if ( sd[4] < 0 ):
                                raise ValueError(
                                    self.__class__.__name__ + ':' + 
                                    ' Invalid species id specification. It was expecting a positive integer, but ' + 
                                    str(sd[4])+ ' was given.'
                                )
                            auxsolutes[0] = sd[4]
                else: 
                    raise ValueError(
                        self.__class__.__name__ + ':' + 
                        ' Invalid specification, exceeds maximum parameters length.'
                    )

            # If a list of aux var specs
            elif isinstance( auxspec, (list,tuple)):
                # List to store aux specs
                auxnames     = []
                # Initialize source variables with default values
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
                                    self.__class__.__name__ + ':' + 
                                    ' Invalid auxiliary variable specification. It was expecting str. ' + 
                                    str(type(aspc))+ ' was given.'
                                )
                            # Good, save it to the list of names
                            auxnames.append( aspc.upper() )
                        elif ( (auxspeclen >= 2 ) and (auxspeclen <= 4 ) ):
                            # Source specification as [auxvar,mass]
                            # Verify str for auxliary variable name 
                            if not isinstance( aspc[0], str ): 
                                raise TypeError(
                                    self.__class__.__name__ + ':' + 
                                    ' Invalid auxiliary variable specification. It was expecting str. ' + 
                                    str(type(aspc[0]))+ ' was given.'
                                )
                            # Good, save it to the list of names
                            auxnames.append( aspc[0].upper() )

                            # Verify particles mass
                            if not isinstance( aspc[1], (int,float) ):
                                raise TypeError(
                                    self.__class__.__name__ + ':' + 
                                    ' Invalid particle mass. It was expecting int/float, but ' + 
                                    str(type(aspc[1]))+ ' was given.'
                                )
                            if ( aspc[1] <= 0 ): 
                                raise ValueError(
                                    self.__class__.__name__ + ':' + 
                                    ' Invalid particle mass. It should be greater than zero. ' +  
                                    str(aspc[1])+ ' was given.'
                                )
                            # Good
                            auxmasses[iaspc] = aspc[1]
                           
                            # Continue reading 
                            if auxspeclen >= 3:
                                # Source specification as [auxvar,mass,template]
                                if isinstance(aspc[2],(list,tuple)):
                                    # Save template
                                    for inp, npa in enumerate(aspc[2]):
                                        auxtemplates[iaspc,inp] = npa
                                elif isinstance( aspc[2], int ):
                                    # Assume uniform template 
                                    auxtemplates[iaspc,:] = aspc[2]
                                else:
                                    raise TypeError(
                                        self.__class__.__name__ + ':' + 
                                        ' Invalid particles template. It was expecting int/list/tuple, but ' + 
                                        str(type(aspc[2]))+ ' was given.'
                                    )
                                # Source specification as [auxvar,mass,template,soluteid]
                                if auxspeclen == 4:
                                    if not isinstance(aspc[3],int):
                                        raise TypeError(
                                            self.__class__.__name__ + ':' + 
                                            ' Invalid species id specification. It was expecting an integer, but ' + 
                                            str(type(aspc[3]))+ ' was given.'
                                        )
                                    if ( aspc[3] < 0 ):
                                        raise ValueError(
                                            self.__class__.__name__ + ':' + 
                                            ' Invalid species id specification. It was expecting a positive integer, but ' + 
                                            str(aspc[3])+ ' was given.'
                                        )
                                    auxsolutes[iaspc] = aspc[3]
                        else:
                            raise ValueError(
                                self.__class__.__name__ + ':' + 
                                ' Invalid specification, maximum length is 4, but ' + 
                                str(len(aspc))+ ' was given.'
                            )
                    elif isinstance( aspc, str ):
                        # Is a str, save it to the list
                        auxnames.append( aspc.upper() )
                    else:
                        raise TypeError(
                            self.__class__.__name__ + ':' + 
                            ' Invalid auxiliary variable specification. It was expecting str. ' + 
                            str(type(aspc))+ ' was given.'
                        )
            else:
                raise TypeError(
                    self.__class__.__name__ + ':' + 
                    ' Invalid auxiliary variable specification. It was expecting type ' + 
                    'str or a list/tuple. '+str(type(auxspec))+ ' was given.'
                )

            # Validate given auxnames
            for ian, an in enumerate(auxnames):
                acount = auxlist.count(an.upper())
                # If the aux name is more than once 
                if acount > 1:
                    raise Exception(
                        self.__class__.__name__ + ':' + 
                        ' Auxiliary name ' + an.upper() + ' was given more than once ' + 
                        'in package ' + pname.upper() + '. Check definition.'
                    )
                # If the given name is not in the list of aux vars of the package
                if (an.upper() not in auxlist):
                    raise Exception(
                        self.__class__.__name__ + ':' + 
                        ' Auxiliary name ' + an.upper() + ' was not found in aux variables ' + 
                        'of package ' + pname.upper() + '. Available variables are: ' + ','.join(auxlist)
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
        ## end read sources ## 

        # Validate that combinations of [sourcename, auxname] are unique 
        counter       =  Counter(frozenset(sl) for sl in sourceslist )
        counterstatus = [counter[frozenset(sl)] == 1 for sl in sourceslist] # True or False

        # Break if some a combination was given more than once
        if not np.all( counterstatus ):
            raise Exception(
                self.__class__.__name__ + ':' + 
                ' Some pairs [package_name,aux_var_name] are repeated. ' +
                ' Combination should be unique for all source packages. '
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


    def _process_spec_format( self, sources=None ):
        '''
        Specific protocols for interpretation of input format based on 
        user specified variables 
        '''
       
        defaultifaceoption = 0
        defaultiface       = 0

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

        # Do some validation for the given sources
        # This format is more rigid, force reading specs from dicts
        if ( not isinstance(sources,(list,dict)) ):
            raise TypeError(
                self.__class__.__name__ + ':' + 
                ' Invalid type for sources in specified input format, expected list or dict. ' + 
                str(type(sources))+ ' was given.'
            )
        if ( isinstance(sources,(dict)) ):
            sources = [sources]
        if ( isinstance(sources,list) ):
            for src in sources:
                if ( not isinstance( src, dict ) ): 
                    raise TypeError(
                        self.__class__.__name__ + ':' + 
                        ' Invalid type for source in specified input format, expected dict. ' + 
                        str(type(src))+ ' was given.'
                    )
        
        # Do more advanced validation for each source specification 
        for src in sources:
            keys = src.keys()

            # Verify mandatory keys
            if ('budgetname' not in keys ):
                raise ValueError(
                    self.__class__.__name__ + ':' + 
                    ' Invalid source specification. The key budgetname was not found in specification.'
                )
            
            # Complete optional keys with defaults
            if ('ifaceoption' not in keys):
                src['ifaceoption' ] = defaultifaceoption 
            if ('defaultiface' not in keys):
                src['defaultiface' ] = defaultiface





        # Save to sources
        self.sources = sources

        # And done
        return
