'''
Tests for the ModpathRWObs class
'''

import pytest
import numpy as np
from flopyrw import modpathrw
from autotest.test_mprw_mt3dp09_cases import MT3DP09Cases

def test_mprw_obs_input_mf6(function_tmpdir):
    '''
    Verifies the input for the observations class
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
   
    # dsp 
    modpathrw.ModpathRWDsp(mp)

    # rwopts
    modpathrw.ModpathRWOpts(mp)

    # simple checks #
    #---------------#
    # by default cellinputoption = 0 (list of cells)
    # by default structured = True (lay,row,col)
    cells = [(0,3,6)]
    with pytest.raises(ValueError):
        # define without list of cells
        modpathrw.ModpathRWObs(mp)
    with pytest.raises(TypeError):
        # define with invalid type of cells
        modpathrw.ModpathRWObs(mp, cells='cells')
    with pytest.raises(ValueError):
        # define with invalid list of structured cells
        modpathrw.ModpathRWObs(mp, cells=[0,1])
    with pytest.raises(ValueError):
        # define with invalid list of structured cells
        modpathrw.ModpathRWObs(mp, cells=[(0,1)])
    with pytest.raises(ValueError):
        # define with unstructured and list of structured cells
        modpathrw.ModpathRWObs(mp, cells=[(0,1,2),(0,1,3)], structured=False)
    with pytest.raises(TypeError):
        # define with unstructured and invalid type list of cells
        modpathrw.ModpathRWObs(mp, cells=['asd',1,2], structured=False)
    with pytest.raises(TypeError):
        # define with invalid type for kind parameter
        modpathrw.ModpathRWObs(mp, kind=None)
    with pytest.raises(ValueError):
        # define with invalid value for kind parameter
        modpathrw.ModpathRWObs(mp, kind=10)
    with pytest.raises(ValueError):
        # define with invalid value for kind parameter
        modpathrw.ModpathRWObs(mp, kind='theobs')

    # obs: define a consistent observation
    modpathrw.ModpathRWObs(mp, cells=cells)
    
    # write distributed # 
    #-------------------#
    mgrid = flowmf6.modelgrid
    cells = np.zeros(shape=(mgrid.nlay,mgrid.nrow,mgrid.ncol), dtype=np.int32)
    cells[0,3,6] = 1
    badcells = np.zeros(shape=(mgrid.nlay,mgrid.nrow,mgrid.ncol), dtype=np.int32)
    badcells[0,3,6] = 1
    badcells[0,3,8] = 2

    with pytest.raises(ValueError):
        # define without the cells array 
        modpathrw.ModpathRWObs(mp, cellinputoption=1)
    with pytest.raises(Exception):
        # define with inconsistent shape for cells array 
        modpathrw.ModpathRWObs(mp, cellinputoption=1, cells=[[0,1],[1,0]])
    with pytest.raises(ValueError):
        # define with inconsistent values for cells array 
        modpathrw.ModpathRWObs(mp, cellinputoption=1, cells=badcells)

    # define consistent for cells array 
    obs = modpathrw.ModpathRWObs(mp, cellinputoption=1, cells=cells)

    # verify assignment to the main model
    pkgs = mp.get_package_list() 
    assert obs._ftype() in pkgs, (
            f"OBS package was not found in ModpathRW object"
        )

    # and write (without checking model consistency)
    mp.write_input()
