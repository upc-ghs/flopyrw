'''
Tests for the ModpathRWIc class
'''

import pytest
import numpy as np
from flopyrw import modpathrw
from autotest.test_mprw_mt3dp09_cases import MT3DP09Cases

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
    with pytest.raises(ValueError):
        # define with invalid kind value
        modpathrw.ModpathRWIc(mp, kind=2, concentration=concentration )
    with pytest.raises(TypeError):
        # define with invalid kind type
        modpathrw.ModpathRWIc(mp, kind=None, concentration=concentration )
    with pytest.raises(ValueError):
        # define with invalid kind string value
        modpathrw.ModpathRWIc(mp, kind='thekind', concentration=concentration )
    with pytest.raises(ValueError):
        # define with invalid particlesdist value
        modpathrw.ModpathRWIc(mp, particlesdist=8, concentration=concentration )
    with pytest.raises(ValueError):
        # define with invalid particlesdist string value
        modpathrw.ModpathRWIc(mp, particlesdist='thedistribution', concentration=concentration )
    with pytest.raises(TypeError):
        # define with invalid particlesdist type
        modpathrw.ModpathRWIc(mp, particlesdist=None, concentration=concentration )


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
    ic = modpathrw.ModpathRWIc(mp2, concentration=concentration )

    # verify assignment to the main model
    pkgs = mp2.get_package_list() 
    assert ic._ftype() in pkgs, (
            f"IC package was not found in ModpathRW object"
        )

    # write (without checks, check=False)
    mp2.write_input()

