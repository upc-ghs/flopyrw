'''
Tests for the ModpathRWBas class
'''

import pytest
import numpy as np
from flopyrw import modpathrw
from autotest.test_mprw_p09mt3d_cases import MT3DP09Cases

def test_mprw_bas_input_mf6(function_tmpdir):
    '''
    Verifies the input for the basic class 
    '''

    # get the mf6 case
    # brings WEL-1 and CHD-1 with aux CONCENTRATION
    flowmf6 = MT3DP09Cases.mf6(function_tmpdir) 

    # modpath-rw
    mp = modpathrw.ModpathRW(
            modelname='mprwsim',
            flowmodel=flowmf6,
            model_ws =function_tmpdir,
        )
   
    # write distributed # 
    #-------------------#
    mgrid = flowmf6.modelgrid
    zeroporosity = np.zeros(shape=(mgrid.nlay,mgrid.nrow,mgrid.ncol), dtype=np.float32)
    badporosity = np.zeros(shape=(mgrid.nlay,mgrid.nrow,mgrid.ncol), dtype=np.float32)
    badporosity[0,3,6] = 2.0 
    badporosity[0,3,8] = -2.0 
    porosity = 0.3*np.ones(shape=(mgrid.nlay,mgrid.nrow,mgrid.ncol), dtype=np.float32)


    with pytest.raises(ValueError):
        # define with only zeros
        modpathrw.ModpathRWBas(mp, porosity=zeroporosity )
    with pytest.raises(ValueError):
        # define with array containing values outside the range
        modpathrw.ModpathRWBas(mp, porosity=badporosity )
    with pytest.raises(ValueError):
        # define with array containing values outside the range
        modpathrw.ModpathRWBas(mp, porosity=-1.0 )
    with pytest.raises(TypeError):
        # define with wrong type
        modpathrw.ModpathRWBas(mp, porosity='porosity' )

    # define consistent package 
    bas = modpathrw.ModpathRWBas(mp, porosity=porosity)

    # verify assignment to the main model
    pkgs = mp.get_package_list()
    assert bas._ftype() in pkgs, (
            f"MPRWBAS package was not found in ModpathRW object"
        )

    # and write (without checking model consistency, check=False by default)
    mp.write_input()
