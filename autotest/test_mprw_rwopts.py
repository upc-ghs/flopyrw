'''
Tests for the ModpathRWOpts class
'''

import pytest
import numpy as np
from flopyrw import modpathrw
from autotest.test_mprw_mt3dp09_cases import MT3DP09Cases

def test_mprw_rwopts_input_mf6(function_tmpdir):
    '''
    Verifies the input for the random walk options class 
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
  
    # timestep
    # define with invalid timestep selection
    with pytest.raises(TypeError):
        modpathrw.ModpathRWOpts(mp, timestep=1)
    with pytest.raises(ValueError):
        modpathrw.ModpathRWOpts(mp, timestep='asd')
    with pytest.raises(TypeError):
        modpathrw.ModpathRWOpts(mp, timestep=None)

    # courant
    # define with invalid courant 
    with pytest.raises(ValueError):
        modpathrw.ModpathRWOpts(mp, courant=0.0)
    with pytest.raises(TypeError):
        modpathrw.ModpathRWOpts(mp, courant='asd')

    # ctdisp
    # define with invalid ctdisp
    with pytest.raises(ValueError):
        modpathrw.ModpathRWOpts(mp, ctdisp=0.0)
    with pytest.raises(TypeError):
        modpathrw.ModpathRWOpts(mp, ctdisp='asd')

    # deltat
    # define with invalid deltat
    with pytest.raises(ValueError):
        modpathrw.ModpathRWOpts(mp, deltat=0.0)
    with pytest.raises(TypeError):
        modpathrw.ModpathRWOpts(mp, deltat='asd')

    # dimensionsmask
    # define with invalid dimensionsmask
    with pytest.raises(ValueError):
        modpathrw.ModpathRWOpts(mp, dimensionsmask=[1,1])
    with pytest.raises(ValueError):
        modpathrw.ModpathRWOpts(mp, dimensionsmask=[1,1,2])
    with pytest.raises(TypeError):
        modpathrw.ModpathRWOpts(mp, dimensionsmask=[1,1,'asd'])

    # advection
    # define with invalid advection selection
    with pytest.raises(TypeError):
        modpathrw.ModpathRWOpts(mp, advection=1)
    with pytest.raises(ValueError):
        modpathrw.ModpathRWOpts(mp, advection='asd')
    with pytest.raises(TypeError):
        modpathrw.ModpathRWOpts(mp, advection=None)

    # randomgenerator
    # define with invalid randomgenerator
    with pytest.raises(TypeError):
        modpathrw.ModpathRWOpts(mp, randomgenerator='asd')
    with pytest.raises(ValueError):
        modpathrw.ModpathRWOpts(mp, randomgenerator=4)

    # define consistent package 
    rwopts = modpathrw.ModpathRWOpts(mp)

    # verify assignment to the main model
    pkgs = mp.get_package_list()
    assert rwopts._ftype() in pkgs, (
            f"RWOPTS package was not found in ModpathRW object"
        )

    # and write (without checking model consistency, check=False by default)
    mp.write_input()
