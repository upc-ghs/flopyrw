'''
Tests for the ModpathRWSpc class
'''

import pytest

from autotest.test_mprw_p09mt3d_cases import MT3DP09Cases
from flopy import modpathrw


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

    # spc: without pgroups
    spc = modpathrw.ModpathRWSpc(mp,dsp)

    with pytest.raises(ValueError):
        # tries to write species without pgroups and particlesmassoption = 0
        # this should be related to check ?
        mp.write_input()

    del spc
    modpathrw.ModpathRWSpc.__class__.INSTANCES = []
    modpathrw.ModpathRWSpc.__class__.COUNTER = 0

    with pytest.raises(TypeError):
        # define a bad spc pgroup param
        spc = modpathrw.ModpathRWSpc(mp,dsp,pgroups='group')

    # define some pgroup indexes
    spc = modpathrw.ModpathRWSpc(mp,dsp,pgroups=[0])

    # and write
    mp.write_input()


