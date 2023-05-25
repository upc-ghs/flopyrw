'''
Tests for the ModpathRWSim class
'''

import pytest
from flopyrw import modpathrw
from autotest.test_mprw_mt3dp09_cases import MT3DP09Cases


def test_mprw_sim_input_mf6(function_tmpdir):
    '''
    Verifies the input for the simulation package 
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
  


    # sim 
    # without: 
    #   - timeseriesoutputoption
    #   - particlesmassoption 
    #   - speciesdispersionoption 
    simconfig = {
        'simulationtype'         : 'rwtimeseries', 
        'trackingdirection'      : 'forward',
        'weaksinkoption'         : 'stop_at',
        'weaksourceoption'       : 'pass_through',
        'referencetime'          : 0.0,
        'stoptimeoption'         : 'specified',
        'stoptime'               : 2*365*86400,
        'timepointdata'          : [2, (1.0*365*86400)],
    }
    # simple checks #
    #---------------#
    with pytest.raises(TypeError):
        # invalid type for timeseriesoutputoption 
        modpathrw.ModpathRWSim(
            mp, 
            **simconfig, 
            timeseriesoutputoption=None,
        )
    with pytest.raises(TypeError):
        # invalid type for timeseriesoutputoption 
        modpathrw.ModpathRWSim(
            mp, 
            **simconfig, 
            timeseriesoutputoption=[0],
        )
    with pytest.raises(ValueError):
        # invalid int value for timeseriesoutputoption 
        modpathrw.ModpathRWSim(
            mp, 
            **simconfig, 
            timeseriesoutputoption=10,
        )
    with pytest.raises(ValueError):
        # invalid str value for timeseriesoutputoption 
        modpathrw.ModpathRWSim(
            mp, 
            **simconfig, 
            timeseriesoutputoption='theoption',
        )
    with pytest.raises(TypeError):
        # invalid type for particlesmassoption 
        modpathrw.ModpathRWSim(
            mp, 
            **simconfig, 
            particlesmassoption=None,
        )
    with pytest.raises(TypeError):
        # invalid type for particlesmassoption 
        modpathrw.ModpathRWSim(
            mp, 
            **simconfig, 
            particlesmassoption=[0],
        )
    with pytest.raises(ValueError):
        # invalid int value for particlesmassoption 
        modpathrw.ModpathRWSim(
            mp, 
            **simconfig, 
            particlesmassoption=10,
        )
    with pytest.raises(ValueError):
        # invalid str value for particlesmassoption 
        modpathrw.ModpathRWSim(
            mp, 
            **simconfig, 
            particlesmassoption='theoption',
        )
    with pytest.raises(TypeError):
        # invalid type for speciesdispersionoption 
        modpathrw.ModpathRWSim(
            mp, 
            **simconfig, 
            speciesdispersionoption=None,
        )
    with pytest.raises(TypeError):
        # invalid type for speciesdispersionoption 
        modpathrw.ModpathRWSim(
            mp, 
            **simconfig, 
            speciesdispersionoption=[0],
        )
    with pytest.raises(ValueError):
        # invalid int value for speciesdispersionoption 
        modpathrw.ModpathRWSim(
            mp, 
            **simconfig, 
            speciesdispersionoption=10,
        )
    with pytest.raises(ValueError):
        # invalid str value for speciesdispersionoption 
        modpathrw.ModpathRWSim(
            mp, 
            **simconfig, 
            speciesdispersionoption='theoption',
        )

    # valid sim 
    simconfig['timeseriesoutputoption']  = 2
    simconfig['particlesmassoption']     = 0
    simconfig['speciesdispersionoption'] = 0
    mprwsim = modpathrw.ModpathRWSim(
        mp, 
        **simconfig
    )
    assert mprwsim.timeseriesoutputoption == 2, (
            f"sim: inconsistent timeseriesoutputoption"
        )
    assert mprwsim.particlesmassoption == 0, (
            f"sim: inconsistent particlesmassoption"
        )
    assert mprwsim.speciesdispersionoption == 0, (
            f"sim: inconsistent speciesdispersionoption"
        )
    simconfig['timeseriesoutputoption']  = 'skip'
    simconfig['particlesmassoption']     = 'off'
    simconfig['speciesdispersionoption'] = 'unique'
    mprwsim = modpathrw.ModpathRWSim(
        mp, 
        **simconfig
    )
    assert mprwsim.timeseriesoutputoption == 2, (
            f"sim: inconsistent timeseriesoutputoption"
        )
    assert mprwsim.particlesmassoption == 0, (
            f"sim: inconsistent particlesmassoption"
        )
    assert mprwsim.speciesdispersionoption == 0, (
            f"sim: inconsistent speciesdispersionoption"
        )


    pkgs = mp.get_package_list() 
    assert mprwsim._ftype() in pkgs, (
            f"MPRWSIM package was not found in ModpathRW object"
        )

    with pytest.raises(Exception):
        # Try to write without rwopts, dsp and bas package
        mp.write_input(check=True)

    # dsp 
    modpathrw.ModpathRWDsp(
        mp,
        alphal=MT3DP09Cases.alphal,
        alphat=MT3DP09Cases.alphat,
        dmeff=MT3DP09Cases.dmeff,
    )

    # bas 
    modpathrw.ModpathRWBas(mp,porosity=MT3DP09Cases.porosity)

    # rwopts
    modpathrw.ModpathRWOpts(mp,dimensionsmask=[1,1,0])
    
    # Try to write ( checking consistency ).
    # consistency should verify if any pgroup ?
    mp.write_input(check=True)
