'''
General utilities 
'''


from .mprw import ModpathRW
def multipackage(initmethod):
    '''
    Wrapper for init methods of packages
    allowing multiple specifications. 

    Manages the storage of multiple instaces
    of the same package on the parent object
    '''

    def wrapper(self, *args, **kwargs):
        # Check model 
        model = args[0]
        if ( not isinstance( model, ModpathRW ) ):
            raise TypeError(
                self.__class__.__name__ + ':' + 
                ' Invalid parent model class. It should be ModpathRW, ' + 
                str(type(model)) + ' was given.'
            )
        # Create the multipackage attribute if not present 
        if ( not hasattr( model, 'multipackage' ) ):
            model.multipackage = {}
        # Get ftype and initialize multipackage dict
        ftype = self._ftype()
        try:
            model.multipackage[ftype]
        except KeyError:
            model.multipackage[ftype] = {
                    'unitnumber': model.next_unit(),
                    'count'     : 0,
                    'instances' : [],
                }

        # Give an id and call initmethod
        # Makes the id consistent with fortran , one-based
        self.id = model.multipackage[ftype]['count'] + 1 
        met     = initmethod(self, *args, **kwargs)

        # If call was successful, then add the package to the parent 
        # and increment counter. If there was an exception, the instance
        # is not stored and everything remains consistent.

        # The difference here w.r.t mf6, is that although there are
        # "multiple" package instances, for modpath-rw is still only
        # one input package definition, containing multiple specifications.
        # So it needs to add the package to the parent only the first
        # time and then store the individual instances somewhere.

        # Add package to the parent only the first time
        if self.parent.multipackage[ftype]['count'] == 0:
            self.parent.add_package(self)
        # Save instance
        self.parent.multipackage[ftype]['instances'].append( self )
        self.parent.multipackage[ftype]['count'] += 1

        return met

    return wrapper
