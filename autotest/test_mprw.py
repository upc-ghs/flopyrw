'''
Tests for the ModpathRW class
'''

import pytest
from modflow_devtools.markers import requires_exe

from autotest.test_mprw_p09mt3d_cases import MT3DP09Cases
from flopy import modpathrw


def test_mprw_nopkgs_mf6(function_tmpdir):
    '''
    Try to write a modpatrw model without mandatory packages
    '''

    # get the mf6 case
    flowmf6 = MT3DP09Cases.mf6(function_tmpdir) 

    # modpath-rw
    mp = modpathrw.ModpathRW(
            modelname='mprwsim',
            flowmodel=flowmf6,
            model_ws =function_tmpdir,
        )
    
    # try to write incomplete #
    #--------------------------
    with pytest.raises(Exception):
        # write without sim, bas, dsp, rwopts
        mp.write_input()
   
    # bas
    modpathrw.ModpathRWBas(
            mp, 
            porosity=0.3
        )
    with pytest.raises(Exception):
        # write without sim, dsp, rwopts
        mp.write_input()

    # sim
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
    with pytest.raises(Exception):
        # write without dsp, rwopts
        mp.write_input()
 
    # rwopts
    modpathrw.ModpathRWOpts(mp)
    with pytest.raises(Exception):
        # write without dsp
        mp.write_input()
    
    # dsp 
    modpathrw.ModpathRWDsp(mp)

    # write    
    mp.write_input()
