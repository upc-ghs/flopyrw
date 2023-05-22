'''
Tests for the ModpathRWSpc class
'''

import numpy as np
import pytest

from autotest.test_mprw_p09mt3d_cases import MT3DP09Cases
from flopy import modpathrw


def test_mprw_imp_input_mf6(function_tmpdir):
    '''
    Verifies the input for the impermeable cells class
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
   
    # simple checks #
    #---------------#
    # defaults:
    #   impformat       = 0
    #   defaultboundary = 0
    #   icbound         = 0
    with pytest.raises(ValueError):
        # define with invalid defaultboundary (0,1) 
        modpathrw.ModpathRWImp(mp, defaultboundary=2)
    with pytest.raises(ValueError):
        # define with invalid impformat (0,1,2) 
        modpathrw.ModpathRWImp(mp, impformat=4)
    with pytest.raises(ValueError):
        # define with impformat=2 and none icbound cells 
        modpathrw.ModpathRWImp(mp, impformat=2, icbound=None)

    # write distributed # 
    #-------------------#
    mgrid = flowmf6.modelgrid
    cells = np.zeros(shape=(mgrid.nlay,mgrid.nrow,mgrid.ncol), dtype=np.int32)
    cells[0,3,6] = 1
    badcells = np.zeros(shape=(mgrid.nlay,mgrid.nrow,mgrid.ncol), dtype=np.int32)
    badcells[0,3,6] = 1
    badcells[0,3,8] = 2

    with pytest.raises(Exception):
        # define with inconsistent shape for cells array 
        modpathrw.ModpathRWImp(mp, impformat=2, cells=[[0,1],[1,0]])
    with pytest.raises(ValueError):
        # define with inconsistent values for icbound 
        modpathrw.ModpathRWImp(mp, impformat=2, icbound=badcells)

    # define consistent package 
    imp  = modpathrw.ModpathRWImp(mp, impformat=2, icbound=cells)

    # verify assignment to the main model
    pkgs = mp.get_package_list() 
    assert imp._ftype() in pkgs, (
            f"IMP package was not found in ModpathRW object"
        )

    # and write (without checking model consistency, check=False by default)
    mp.write_input()
