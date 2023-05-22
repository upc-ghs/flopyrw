'''
Tests for the ModpathRWSpc class
'''

import numpy as np
import pytest

from autotest.test_mprw_p09mt3d_cases import MT3DP09Cases
from flopy import modpathrw


def test_mprw_ic_input_mf6(function_tmpdir):
    '''
    Verifies the input for the initial conditions class
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
    with pytest.raises(ValueError):
        # define with None concentration
        modpathrw.ModpathRWIc(mp, concentration=None)
    with pytest.raises(TypeError):
        # define with invalid concentration list 
        modpathrw.ModpathRWIc(mp, concentration=['asd'])
    with pytest.raises(Exception):
        # define with invalid concentration shape 
        modpathrw.ModpathRWIc(mp, concentration=[2,1])

    # create a valid concentration
    mgrid = flowmf6.modelgrid
    concentration = np.ones(shape=(mgrid.nlay,mgrid.nrow,mgrid.ncol))

    with pytest.raises(ValueError):
        # define with zero mass 
        modpathrw.ModpathRWIc(mp, concentration=concentration, particlesmass=0.0)
    with pytest.raises(ValueError):
        # define with negative mass 
        modpathrw.ModpathRWIc(mp, concentration=concentration, particlesmass=-10.0)
    with pytest.raises(ValueError):
        # define with invalid species id 
        modpathrw.ModpathRWIc(mp, concentration=concentration, speciesid=-9)
    with pytest.raises(TypeError):
        # define with invalid species id 
        modpathrw.ModpathRWIc(mp, concentration=concentration, speciesid='aspecies' )

    # sim with particlesmassoption = 2
    simconfig = {
        'simulationtype'         : 'rwtimeseries', 
        'trackingdirection'      : 'forward',
        'weaksinkoption'         : 'stop_at',
        'weaksourceoption'       : 'pass_through',
        'referencetime'          : 0.0,
        'stoptimeoption'         : 'specified',
        'stoptime'               : 2*365*86400,
        'particlesmassoption'    : 2,
        'speciesdispersionoption': 0,
        'timeseriesoutputoption' : 2, 
        'timepointdata'          : [2, (1.0*365*86400)],
    }
    mpsim = modpathrw.ModpathRWSim(
        mp,
        **simconfig
    )
    # ic: by default speciesid = 0
    modpathrw.ModpathRWIc(mp, concentration=concentration )

    # write
    mp.write_input()


    # model with particlesmassoption = 0 #
    #------------------------------------#
    # modpath-rw
    mp2 = modpathrw.ModpathRW(
            modelname='mprwsim_2',
            flowmodel=flowmf6,
            model_ws =function_tmpdir,
        )
   
    # dsp 
    modpathrw.ModpathRWDsp(mp2)

    # rwopts
    modpathrw.ModpathRWOpts(mp2)

    simconfig['particlesmassoption'] = 0
    mpsim = modpathrw.ModpathRWSim(
        mp2,
        **simconfig
    )

    # ic
    modpathrw.ModpathRWIc(mp2, concentration=concentration )

    # write
    mp2.write_input()

