'''
Tests for the ModpathRWSpc class
'''

import pytest
from flopyrw import modpathrw
from autotest.test_mprw_mt3dp09_cases import MT3DP09Cases

def test_mprw_spc_input_mf6(function_tmpdir):
    '''
    Verifies the input for species class
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
    dsp = modpathrw.ModpathRWDsp(mp)

    # rwopts
    modpathrw.ModpathRWOpts(mp)

    # simple checks #
    #---------------#
    with pytest.raises(TypeError):
        # define without a dispersion model
        modpathrw.ModpathRWSpc(mp)
    with pytest.raises(TypeError):
        # define with an invalid dispersion model
        modpathrw.ModpathRWSpc(mp,'asd')


    # sim with particlesmassoption = 0
    simconfig = {
        'simulationtype'         : 'rwtimeseries', 
        'trackingdirection'      : 'forward',
        'weaksinkoption'         : 'stop_at',
        'weaksourceoption'       : 'pass_through',
        'referencetime'          : 0.0,
        'stoptimeoption'         : 'specified',
        'stoptime'               : 2*365*86400,
        'particlesmassoption'    : 0, 
        'speciesdispersionoption': 0,
        'timeseriesoutputoption' : 2, 
        'timepointdata'          : [2, (1.0*365*86400)],
    }
    mpsim = modpathrw.ModpathRWSim(
        mp,
        **simconfig
    )
    
    # spc: without pgroups and particlesmassoption = 0
    spc = modpathrw.ModpathRWSpc(mp,dsp)
    with pytest.raises(ValueError):
        # this should be related to check ?
        mp.write_input()

    # The multipackage variables are somehow 
    # inconsistent now. Due to the exception 
    # thrown at the write method. However, 
    # this occurs only in the context of testing.

    # Force reset multipackage
    mp.multipackage =  {}
    # and add dsp again
    dsp = modpathrw.ModpathRWDsp(mp)

    with pytest.raises(TypeError):
        # define a bad spc pgroup param
        spc = modpathrw.ModpathRWSpc(mp,dsp,pgroups='group')

    with pytest.raises(TypeError):
        # define a bad spc pgroup param
        spc = modpathrw.ModpathRWSpc(mp,dsp,pgroups=['group'])

    # define some pgroup indexes
    spc = modpathrw.ModpathRWSpc(mp,dsp,pgroups=[0])

    # and write
    mp.write_input()


    # with particlesmassoption!=0,1 #
    # and speciesdispersionoption=1 #
    #-------------------------------#
    # modpath-rw
    mp2 = modpathrw.ModpathRW(
            modelname='mprwsim_2',
            flowmodel=flowmf6,
            model_ws =function_tmpdir,
        )
   
    # dsp 
    dsp = modpathrw.ModpathRWDsp(mp2)

    # rwopts
    modpathrw.ModpathRWOpts(mp2)
    
    # sim 
    simconfig['particlesmassoption'] = 2
    simconfig['speciesdispersionoption'] = 1
    modpathrw.ModpathRWSim(
        mp2,
        **simconfig
    )
    
    # spc
    spc = modpathrw.ModpathRWSpc(mp2,dsp)

    # verify assignment to the main model
    pkgs = mp2.get_package_list() 
    assert spc._ftype() in pkgs, (
            f"SPC package was not found in ModpathRW object"
        )

    # Write
    mp2.write_input()

