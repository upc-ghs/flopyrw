'''
Tests for the ModpathRWSim class
'''

import pytest
from flopyrw import modpathrw
from autotest.test_mprw_mt3dp09_cases import MT3DP09Cases
from modflow_devtools.markers import requires_exe 


@requires_exe("mf6","mpathrw")
def test_mprw_sim_run_mf6(function_tmpdir):
    '''
    Verifies running of the simulation
    '''

    # get the mf6 case
    # brings WEL-1 and CHD-1 with aux CONCENTRATION
    flowmf6 = MT3DP09Cases.mf6(function_tmpdir,write=True,run=True)

    # modpath-rw
    mp = modpathrw.ModpathRW(
            modelname='mprwsim',
            flowmodel=flowmf6,
            model_ws =function_tmpdir,
        )
 
    # bas 
    modpathrw.ModpathRWBas(mp,porosity=MT3DP09Cases.porosity)

    # rwopts
    modpathrw.ModpathRWOpts(
        mp,
        dimensionsmask=[1,1,0],
        timestep='min',
        ctdisp=0.1,
        courant=0.1
    )

    # src
    sources = [
        (
            "WEL-1",
            [
                ["CONCENTRATION", 200.0, (4,4,1)],
            ],
        ),
    ]
    src = modpathrw.ModpathRWSrc(
        mp,
        inputformat='aux', 
        sources=sources,
    )

    # dsp 
    modpathrw.ModpathRWDsp(
        mp,
        alphal=MT3DP09Cases.alphal,
        alphat=MT3DP09Cases.alphat,
        dmeff=MT3DP09Cases.dmeff,
    )


    # obs
    obs = modpathrw.ModpathRWObs(
        mp,
        kind=1, # flux
        cells=[tuple(MT3DP09Cases.extwell)],
    )

    # sim
    simconfig = {
        'simulationtype'         : 'rwtimeseries', 
        'trackingdirection'      : 'forward',
        'weaksinkoption'         : 'stop_at',
        'weaksourceoption'       : 'pass_through',
        'referencetime'          : 0.0,
        'stoptimeoption'         : 'specified',
        'stoptime'               : 2*365*86400,
        'timepointdata'          : [73, (1.0*10*86400)],
        'timeseriesoutputoption' : 2,
        'particlesmassoption'    : 0,
        'speciesdispersionoption': 0,
    }
    mprwsim = modpathrw.ModpathRWSim(
        mp, 
        **simconfig
    )


    # Try to write ( checking consistency ).
    # consistency should verify if any pgroup ?
    mp.write_input(check=True)

    success, buff = mp.run_model(silent=True,report=True)
    assert success, f"mpathrw did not run correctly"
