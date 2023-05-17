'''
Configuration of MODPATH-RW basic package
'''

# flopy 
from flopy.modpath import Modpath7Bas


class ModpathRWBas( Modpath7Bas ):
    '''
    MODPATH-RW Basic class 

    Extends flopy.modpath.Modpath7Bas
    '''

    def __init__(self, *args, **kwargs):
        # Call parent constructor
        super().__init__(*args,**kwargs)

        # Not the most elegant solution but consistent
        self.name  = ['MPRWBAS']
        self._name = ['MPRWBAS']
        self._generate_heading()
