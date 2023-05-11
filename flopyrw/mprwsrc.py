'''
Configuration of MODPATH-RW sources
'''

# python
from collections import Counter
import numpy as np

# flopy
from flopy.pakbase import Package
from flopy.utils import Util3d

# local 
from .utils import count_instances # Increment COUNTER


class ModpathRWSrc( Package ):
    """
    MODPATH-RW Source Package Class.

    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.modpath.Modpath7`) to which
        this package will be added.
    """


    # Class properties
    COUNTER = 0
    UNITNUMBER = 0
    INSTANCES = []


    @count_instances
    def __init__(
        self,
        model,
        format        = 'aux',
        sources       = None,
        cells         = None, 
        mass          = 1.0,
        nparticles    = 2,
        concentration = None, 
        soluteid      = None,
        sourcename    = None, 
        auxiliary     = None,
        id            = None, 
        stringid      = None, 
        extension     = 'src',
    ):
       
        # Define UNITNUMBER if the first instance is created
        if self.__class__.COUNTER == 1:
            unitnumber = model.next_unit()
            super().__init__(model, extension, "SRC", unitnumber)
            self.__class__.UNITNUMBER = self.unit_number[0]
        else:
            # Needed for Util3d
            self._parent = self.INSTANCES[0]._parent

        shape = model.shape
        if len(shape) == 3:
            shape3d = shape
        elif len(shape) == 2:
            shape3d = (shape[0], 1, shape[1])
        else:
            shape3d = (1, 1, shape[0])
        self.model = model
        self.shape3d = shape3d

     
        self.format = format
        self.modelversion = self._parent.flowmodel.version 
       

        if mass is not None:
            DEFAULTMASS = mass
        else: 
            DEFAULTMASS = 1.0

        if nparticles is not None:
            DEFAULTNP = nparticles
        else:
            DEFAULTNP = 2

        # Clarify whether this is python based or fortran based
        if soluteid is not None:
            DEFAULTSOLID = int(soluteid)
        else:
            DEFAULTSOLID = int(1)


        # If format is aux, validate sources, verify auxiliary variables 
        sourceslist = [] # STORE IN CLASS ?
        alldatalist = []

        # Validate format kind
        if ( ( format.upper() not in ['AUX','AUXILIARY'] ) ):
            raise Exception('flopyrw:ModpathRWSrc: source format ' + format.upper() + ' not implemented. ')

        # Validate that sources was given
        if ( ( sources is None ) and ( format.upper() in ['AUX','AUXILIARY'] ) ):
            raise Exception('flopyrw:ModpathRWSrc: package requires sources data. None was given.')

        # Read sources
        for sd in sources:

            pname = sd[0]
            if not isinstance( pname, str ): 
                raise TypeError('flopyrw:ModpathRWSrc: package name should be str. ',type(pname),' was given.')
            pkg = self._parent.flowmodel.get_package(pname.lower())
            if pkg is None:
                raise Exception(\
                    'flopyrw:ModpathRWSrc: package name '+ pname.lower() \
                    + ' was not found in flowmodel.'+' Available packages are: ' +",".join(self._parent.flowmodel.get_package_list()) )

            # Get the auxlist from package
            # Note: auxlist filled only with str.upper()

            if self.modelversion == 'mf6':
                # Get the aux data and create auxlist
                if not pkg.auxiliary.has_data():
                    raise Exception('flopyrw:ModpathRWSrc: auxiliary property does not has data for package ' + pkg.name[0])
                aux     = pkg.auxiliary.array[0]
                auxlist = [a.upper() for a in aux ] 
            else:
                # Extract possible aux variable names and perform some preliminary checks for != MF6
                opts = pkg.options
                if isinstance(opts, flopy.utils.OptionBlock):
                    auxlistcand = opts.auxillary
                elif isinstance(opts,(list,tuple)):
                    auxlistcand = pkg.options
                else:
                    raise TypeError('flopyrw:ModpathRWSrc: unexpected options type in package ' + pkg.name[0] )

                if ( len(auxlistcand) == 0):
                    raise Exception('flopyrw:ModpathRWSrc: no auxiliary variables were found in package ' + pkg.name[0] )

                # Loop over possible aux vars and create auxlist
                auxlist = []
                for ia, auxc in enumerate(auxlistcand):
                    # Clean the name ( str is assumed )

                    # Clean borders
                    cand = auxc.strip().upper()
                    
                    # If is the identifier str, go the next
                    if ( cand in ['AUX','AUXILIARY'] ):
                        # To the next
                        continue

                    # Break it
                    cand = cand.split()

                    # It it had a whitespace
                    if ( len(cand) == 2 ):
                        # Need to verify that the first word 
                        # is the aux identifier and the second
                        # is the auxvar name
                        if ( cand[0].upper() in ['AUX','AUXILIARY'] ):
                            # Alles gud
                            auxlist.append( cand[1].upper() )
                    elif ( len(cand) == 1):
                        auxlist.append( cand[0].upper() )
                    else:
                        # Probably the case in which auxvars have whitespace in between. omg
                        print( cand )
                        raise Exception('flopyrw:ModpathRWSrc: something wrong with definition of aux in package ' + pkg.name[0])

            # Determine the format in which the aux variables are specified
            auxspec = sd[1]
            if isinstance( auxspec, str ):
                # If a str was given, transform to list
                auxnames     = [auxspec.upper()]
                auxmasses    = [DEFAULTMASS]
                auxtemplates = DEFAULTNP*np.ones(shape=(1,3),dtype=np.int32)
                auxsolutes   = [DEFAULTSOLID] # Clarify whether these are python based or fortran based
            elif isinstance( auxspec, (list,tuple)):
                # List to store aux specs
                auxnames     = []
                auxmasses    = DEFAULTMASS*np.ones( shape=(len(auxspec),) )
                auxtemplates = DEFAULTNP*np.ones(shape=(len(auxspec),3),dtype=np.int32)
                auxsolutes   = DEFAULTSOLID*np.ones( shape=(len(auxspec),),dtype=np.int32 ) # Clarify whether these are python based or fortran based
                for iaspc, aspc in enumerate(auxspec):
                    if isinstance( aspc, (list,tuple) ):
                        auxspeclen = len(aspc)
                        if auxspeclen == 1:
                            # Consider as auxvar
                            if not isinstance( aspc, str ): 
                                raise TypeError('flopyrw:ModpathRWSrc: auxiliary name should be str. ',type(aspc),' was given.')
                            # Good
                            auxnames.append( aspc.upper() )
                        elif ( (auxspeclen >= 2 ) and (auxspeclen <= 4 ) ):
                            # First is auxvarname
                            if not isinstance( aspc[0], str ): 
                                raise TypeError('flopyrw:ModpathRWSrc: auxiliary name should be str. ',type(aspc[0]),' was given.')
                            # Good
                            auxnames.append( aspc[0].upper() )
                            
                            # Second is mass 
                            if not isinstance( aspc[1], (int,float) ): 
                                raise TypeError('flopyrw:ModpathRWSrc: particles mass should be int/float. ',type(aspc[1]),' was given.')
                            # Good
                            auxmasses[iaspc] = aspc[1]
                           
                            # Continue reading 
                            if auxspeclen >= 3:
                                # Consider as [auxvar,mass,template]
                                if isinstance(aspc[2],(list,tuple)):
                                    for inp, npa in enumerate(aspc[2]):
                                        auxtemplates[iaspc,inp] = npa
                                elif isinstance( aspc[2], int ):
                                    # Assume uniform template 
                                    auxtemplates[iaspc,:] = aspc[2]
                                else:
                                    raise TypeError('flopyrw:ModpathRWSrc: particles template should be int/list/tuple. ',type(aspc[2]),' was given.')

                                if auxspeclen == 4:
                                    # Consider as [auxvar,mass,template,soluteid]
                                    if not isinstance(aspc[3],int):
                                        raise TypeError('flopyrw:ModpathRWSrc: solute id should be int. ',type(aspc[3]),' was given.')
                                    auxsolutes[iaspc] = aspc[3]
                        else:
                            raise ValueError('flopyrw:ModpathRWSrc: max length of aux specification is 4. ' + str(len(aspc)) + ' was given.')
                    elif isinstance( aspc, str ):
                        # Is a str, give to the list
                        auxnames.append( aspc.upper() )
                    else: 
                        raise TypeError('flopyrw:ModpathRWSrc: auxiliary name should be str. ',type(aspc),' was given.')

            else:
                raise TypeError('flopyrw:ModpathRWSrc: auxiliary name should be str or list/tuple. ',type(auxspec),' was given.')


            # Validate given auxnames
            for ian, an in enumerate(auxnames):
                acount = auxlist.count(an.upper())
                if acount > 1:
                    raise Exception('flopyrw:ModpathRWSrc: auxiliary name '+an.upper()+' was specified more than once in package ' + pname.upper() +'. Check definition.')
                if an.upper() not in auxlist:
                    raise ValueError('flopyrw:ModpathRWSrc: auxiliary name ' +an.upper()+' was not specified in package '+ pname.upper() +'. Check package definition.') 

                # If it got until here, definition is consistent,
                # save the pair [sourcename,auxname]
                tempair  = [pname.upper(),an.upper()]
                sourceslist.append(tempair)
                tempdata = [an.upper(), auxmasses[ian], auxtemplates[ian,0], auxtemplates[ian,1], auxtemplates[ian,2], auxsolutes[ian] ]
                alldatalist.append(tempdata)
        ## end read sources ## 

        # Validate that combinations of [sourcename, auxname] are unique 
        counter       =  Counter(frozenset(sl) for sl in sourceslist )
        counterstatus = [counter[frozenset(sl)] == 1 for sl in sourceslist]

        # Break if some inconsistencies
        if not np.all( counterstatus ):
            raise Exception('flopyrw:ModpathRWSrc: some pairs of [sourcename,auxname] are repeated. Combination should be unique for all source packages.')

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


        # Create format        
        fmts = []
        fmts.append("{:20s}")  # auxvarname
        fmts.append("{:.6f}")  # particles mass
        for nt in range(3):    # (nx,ny,nz)
            fmts.append("{:3d}")
        # Append solute id ?
        istart = 0
        iend   = 5
        if ( 
            ( self.INSTANCES[0]._parent.particlesmassoption == 2 ) or  # If solute id for particle groups
            ( self.INSTANCES[0]._parent.particlesmassoption == 2 ) ):  # Solute specific dispersion
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

            # Write id's
            f.write(f"{ins.stringid}\n")

            # Write format
            f.write(f"{ins.format.upper()}\n")

            if (
                ( ins.format.upper() == 'AUX' ) or 
                ( ins.format.upper() == 'AUXILIARY' ) ):

                # Inform about the number of source budgets
                f.write(f"{len(ins.uniquesources)}\n")

                # (...) and loop over them
                for isrc, src in enumerate(ins.uniquesources):

                    # Write src name
                    f.write(f"{src.upper()}\n")

                    # Write the number of aux vars
                    f.write(f"{len(ins.allauxnames[isrc])}\n")
                    
                    # Write aux vars
                    for iaux, auxdata in enumerate(ins.alldataperaux[isrc]):
                        f.write(fmt.format(*auxdata[istart:iend]))


        # And close
        f.close()


        # Done
        return
