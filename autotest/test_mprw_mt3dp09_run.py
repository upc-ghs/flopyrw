'''
Tests for the ModpathRWSim class
'''

import pytest
from flopyrw import modpathrw
from autotest.test_mprw_mt3dp09_cases import MT3DP09Cases
from modflow_devtools.markers import requires_exe 


@requires_exe("mf6","mpathrw")
def test_mprw_run_tsobs_mf6(function_tmpdir):
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
                ["CONCENTRATION", 300.0, (4,4,1)],
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



@requires_exe("mf6","mpathrw")
def test_mprw_run_tsgpkde_mf6(function_tmpdir):
    '''
    Verifies running of the simulation with reconstruction
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
                ["CONCENTRATION", 300.0, (4,4,1)],
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

    # gpkde
    modpathrw.ModpathRWGpkde(mp)

    # sim
    simconfig = {
        'simulationtype'         : 'rwtimeseries', 
        'trackingdirection'      : 'forward',
        'weaksinkoption'         : 'stop_at',
        'weaksourceoption'       : 'pass_through',
        'referencetime'          : 0.0,
        'stoptimeoption'         : 'specified',
        'stoptime'               : 2*365*86400,
        'timepointdata'          : [2, (1.0*365*86400)],
        'timeseriesoutputoption' : 0,
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


@requires_exe("mf6","mpathrw")
def test_mprw_run_epoint_mf6(function_tmpdir):
    '''
    Verifies running of endpoint simulation
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
                ["CONCENTRATION", 250.0, (4,4,1)],
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

    # sim
    simconfig = {
        'simulationtype'         : 'rwendpoint', 
        'trackingdirection'      : 'forward',
        'weaksinkoption'         : 'stop_at',
        'weaksourceoption'       : 'pass_through',
        'referencetime'          : 0.0,
        'stoptimeoption'         : 'specified',
        'stoptime'               : 1*365*86400,
        'timeseriesoutputoption' : 0,
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


@requires_exe("mf6","mpathrw")
def test_mprw_run_combined_mf6(function_tmpdir):
    '''
    Verifies running of combined simulation
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
                ["CONCENTRATION", 500.0, (3,3,1)],
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

    # sim
    simconfig = {
        'simulationtype'         : 'rwcombined', 
        'trackingdirection'      : 'forward',
        'weaksinkoption'         : 'stop_at',
        'weaksourceoption'       : 'pass_through',
        'referencetime'          : 0.0,
        'stoptimeoption'         : 'specified',
        'stoptime'               : 1*365*86400,
        'timepointdata'          : [36, (1.0*10*86400)],
        'timeseriesoutputoption' : 0,
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


@requires_exe("mf6","mpathrw")
def test_mprw_run_tsobs_mf6disvusg(function_tmpdir):
    '''
    Verifies running of the simulation
    '''

    # get the mf6 case
    # brings WEL-1 and CHD-1 with aux CONCENTRATION
    # appends usgextnodes and usginjnodes to MT3DP09Cases
    flowmf6 = MT3DP09Cases.mf6disvusg(function_tmpdir,write=True,run=True)

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
                ["CONCENTRATION", 400.0, (1,1,1)],
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
        kind='flux',
        cells=[ no for no in MT3DP09Cases.usgextnodes ],
        structured=False,
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


@requires_exe("mf6","mpathrw")
def test_mprw_run_tsobs_mf6disvmultilayer(function_tmpdir):
    '''
    Verifies running of the simulation

    This test might be slow due to the large dispersivities
    with respect to layer height, and eventually particles 
    modify their time step several times to avoid large relative 
    jumps.
    '''

    # get the mf6 case
    # brings WEL-1 and CHD-1 with aux CONCENTRATION
    # appends linearextnodes and linearinjnodes to MT3DP09Cases
    flowmf6 = MT3DP09Cases.mf6disvmultilayer(function_tmpdir,write=True,run=True)

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
        dimensionsmask=[1,1,1], # 3d
        timestep='min',
        ctdisp=0.1,
        courant=0.1
    )

    # src
    sources = [
        (
            "WEL-1",
            [
                ["CONCENTRATION", 500.0, (2,2,1)],
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
        kind='flux',
        cells=[ no for no in MT3DP09Cases.linearextnodes ],
        structured=False,
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



@requires_exe("mf6","mpathrw")
def test_mprw_run_tsobs_mf6timeseriesaux(function_tmpdir):
    '''
    Verifies running of mpathrw simulation 
    with timeseries aux variable.
    '''

    # get the mf6 case
    # brings WEL-1 and CHD-1 with aux CONCENTRATION
    flowmf6 = MT3DP09Cases.mf6timeseriesaux(function_tmpdir,write=True,run=True)

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
                ["CONCENTRATION", 300.0, (4,4,1)],
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
        kind ='flux',
        cells=[tuple(MT3DP09Cases.extwell)],
        outputoption = 1, 
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


@requires_exe("mf6","mpathrw")
def test_mprw_run_epointextendobs_mf6(function_tmpdir):
    '''
    Verifies running of endpoint simulation
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
        dimensionsmask= [1,1,0],
        timestep      = 'min',
        ctdisp        = 0.1,
        courant       = 0.1
    )

    # src
    sources = [
        (
            "WEL-1",
            [
                ["CONCENTRATION", 300.0, (4,4,1)],
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
        kind='flux',
        cells=[tuple(MT3DP09Cases.extwell)],
    )

    # sim
    simconfig = {
        'simulationtype'         : 'rwendpoint', 
        'trackingdirection'      : 'forward',
        'weaksinkoption'         : 'stop_at',
        'weaksourceoption'       : 'pass_through',
        'referencetime'          : 0.0,
        'stoptimeoption'         : 'extend',
        'endpointoutputoption'   : 1,
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



@requires_exe("mf6","mpathrw")
def test_mprw_run_tsobs_backward_mf6(function_tmpdir):
    '''
    Verifies running of the simulation
    '''

    # get the mf6 case
    # brings WEL-1 and CHD-1 with aux CONCENTRATION
    flowmf6 = MT3DP09Cases.mf6(function_tmpdir,write=True,run=True, backward=True)

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
                ["CONCENTRATION", 300.0, (4,4,1)],
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
        cells=[tuple(MT3DP09Cases.injwell)], # observe the injwell
    )

    # sim
    simconfig = {
        'simulationtype'         : 'rwtimeseries', 
        'trackingdirection'      : 'backward',
        'weaksinkoption'         : 'pass_through',
        'weaksourceoption'       : 'stop_at',
        'referencetime'          : 2.0*365*86400,
        'stoptimeoption'         : 'total',
        'timepointdata'          : [36, (2.0*10*86400)],
        'timeseriesoutputoption' : 0,
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


@requires_exe("mf6","mpathrw")
def test_mprw_run_tsobs_backward_mf6timeseriesaux(function_tmpdir):
    '''
    Verifies running of mpathrw simulation 
    with timeseries aux variable.
    '''

    # get the mf6 case
    # brings WEL-1 and CHD-1 with aux CONCENTRATION
    flowmf6 = MT3DP09Cases.mf6timeseriesaux(function_tmpdir,write=True,run=True, backward=True)

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
                ["CONCENTRATION", 300.0, (4,4,1)],
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
        kind='flux', 
        cells=[tuple(MT3DP09Cases.injwell)], # observe the injwell
    )

    # sim
    simconfig = {
        'simulationtype'         : 'rwtimeseries', 
        'trackingdirection'      : 'backward',
        'weaksinkoption'         : 'pass_through',
        'weaksourceoption'       : 'stop_at',
        'referencetime'          : 2.0*365*86400,
        'stoptimeoption'         : 'total',
        'timepointdata'          : [36, (2.0*10*86400)],
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



@requires_exe("mf6","mpathrw")
def test_mprw_run_srcspec_mf6timeserieslinear(function_tmpdir):
    '''
    Verifies running of mpathrw simulation 
    with timeseries aux variable.
    '''

    # get the mf6 case
    # brings WEL-1 and CHD-1
    flowmf6 = MT3DP09Cases.mf6timeserieslinear(function_tmpdir,write=True,run=True)

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
        {
            'budgetname': 'WEL', 
            'timeintervals' : [
                [0.0,365*86400.0],
            ],
            'cellinput': 0,
            'concentration': [MT3DP09Cases.cinjwell],
            'particlesmass': 300, 
            'template':[(4,4,1)],
            'speciesid': 0,
            'cells': None,
        },
    ]
    src = modpathrw.ModpathRWSrc(
        mp,
        inputformat='spec', 
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
        kind ='flux',
        cells=[tuple(MT3DP09Cases.extwell)],
        outputoption = 1, 
    )
    # sim
    simconfig = {
        'simulationtype'         : 'rwtimeseries', 
        'trackingdirection'      : 'forward',
        'weaksinkoption'         : 'stop_at',
        'weaksourceoption'       : 'pass_through',
        'referencetime'          : 0.0,
        'stoptimeoption'         : 'specified',
        'stoptime'               : 4*86400,
        'timepointdata'          : [36, (2.0*10*86400)],
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


@requires_exe("mf6","mpathrw")
def test_mprw_run_srcspec_mf6(function_tmpdir):
    '''
    Verifies running of mpathrw simulation 
    with source using spec format. 
    '''

    # get the mf6 case
    # brings WEL-1 and CHD-1
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
        {
            'budgetname': 'WEL', 
            'timeintervals' : [
                [0.0,365*86400.0],
            ],
            'cellinput': 0,
            'concentration': [MT3DP09Cases.cinjwell],
            'particlesmass': 300, 
            'template':[(4,4,1)],
            'speciesid': 0,
            'cells': None,
        },
    ]
    src = modpathrw.ModpathRWSrc(
        mp,
        inputformat='spec', 
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
        kind ='flux',
        cells=[tuple(MT3DP09Cases.extwell)],
        outputoption = 1, 
    )
    # sim
    simconfig = {
        'simulationtype'         : 'rwendpoint', 
        'trackingdirection'      : 'forward',
        'weaksinkoption'         : 'stop_at',
        'weaksourceoption'       : 'pass_through',
        'referencetime'          : 0.0,
        'stoptimeoption'         : 'specified',
        'stoptime'               : 2*365*86400.0,
        'endpointoutputoption'   : 2,
    }
    mprwsim = modpathrw.ModpathRWSim(
        mp, 
        **simconfig
    )


    # Try to write ( checking consistency ).
    mp.write_input(check=True)

    success, buff = mp.run_model(silent=True,report=True)
    assert success, f"mpathrw did not run correctly"



@requires_exe("mf6","mpathrw")
def test_mprw_run_srcaux_mf6timeserieslinear_backward(function_tmpdir):
    '''
    Verifies running of mpathrw simulation 
    with timeseries aux variable.
    '''

    # get the mf6 case
    # brings WEL-1 and CHD-1
    flowmf6 = MT3DP09Cases.mf6timeserieslinear(function_tmpdir,write=True,run=True, backward=True, aux=True)

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
                ["CONCENTRATION", 300.0, (4,4,1)],
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
        kind ='flux',
        cells=[tuple(MT3DP09Cases.injwell)],
        outputoption = 1, 
    )
    # sim
    simconfig = {
        'simulationtype'         : 'rwtimeseries', 
        'trackingdirection'      : 'backward',
        'weaksinkoption'         : 'pass_through',
        'weaksourceoption'       : 'stop_at',
        'referencetime'          : 2.0*365*86400,
        'stoptimeoption'         : 'total',
        'timepointdata'          : [36, (2.0*10*86400)],
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


@requires_exe("mf6","mpathrw")
def test_mprw_run_srcspec_mf6timeserieslinear_backward(function_tmpdir):
    '''
    Verifies running of mpathrw simulation 
    with timeseries aux variable.
    '''

    # get the mf6 case
    # brings WEL-1 and CHD-1
    flowmf6 = MT3DP09Cases.mf6timeserieslinear(function_tmpdir,write=True,run=True)

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
        {
            'budgetname': 'WEL', 
            'timeintervals' : [
                [365*86400.0,727*86400.0],
                [728*86400.0,730*86400.0],
            ],
            'cellinput': 0,
            'concentration': [0.035*MT3DP09Cases.cinjwell, 0.025*MT3DP09Cases.cinjwell],
            'particlesmass': 300, 
            'template':[(4,4,1)],
            'speciesid': 0,
            'cells': None,
        },
    ]
    src = modpathrw.ModpathRWSrc(
        mp,
        inputformat='spec', 
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
        kind ='flux',
        cells=[tuple(MT3DP09Cases.injwell)],
        outputoption = 1, 
    )
    # sim
    simconfig = {
        'simulationtype'         : 'rwtimeseries', 
        'trackingdirection'      : 'backward',
        'weaksinkoption'         : 'pass_through',
        'weaksourceoption'       : 'stop_at',
        'referencetime'          : 2.0*365*86400,
        'stoptimeoption'         : 'total',
        'timepointdata'          : [36, (2.0*10*86400)],
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
