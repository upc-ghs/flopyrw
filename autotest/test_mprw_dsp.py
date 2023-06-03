'''
Tests for the ModpathRWDsp class
'''

import pytest
from flopyrw import modpathrw
from autotest.test_mprw_mt3dp09_cases import MT3DP09Cases


def test_mprw_dsp_input_mf6(function_tmpdir):
    '''
    Verifies the input for dispersion package
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

    # simple checks #
    #---------------#
    with pytest.raises(TypeError):
        # pass an int to inputformat
        modpathrw.ModpathRWDsp(
            mp,
            inputformat=0
        )
    with pytest.raises(ValueError):
        # pass an invalid input format 
        modpathrw.ModpathRWDsp(
            mp,
            inputformat='auxasd',
        )
    with pytest.raises(TypeError):
        # pass a None inputformat and None dispersion parameter.
        # first verifies health of dispersion parameters
        modpathrw.ModpathRWDsp(
            mp,
            inputformat=None,
            alphal=None,
        )
    with pytest.raises(TypeError):
        # pass an int to modelkind
        modpathrw.ModpathRWDsp(
            mp,
            modelkind=0, 
        )
    with pytest.raises(ValueError):
        # pass an invalid value for modelkind
        modpathrw.ModpathRWDsp(
            mp,
            modelkind='themostlinear', 
        )
    
    # try to define invalid dispersion params #
    #-----------------------------------------#
    with pytest.raises(Exception):
        # define with inconsistent shape for alphal array 
        modpathrw.ModpathRWDsp(mp, alphal=[[0.1,0.1],[1,0]])
    with pytest.raises(Exception):
        # define with inconsistent shape for alphat array 
        modpathrw.ModpathRWDsp(mp, alphat=[[0.1,0.1],[1,0]])
    with pytest.raises(Exception):
        # define with inconsistent shape for dmeff array 
        modpathrw.ModpathRWDsp(mp, dmeff=[[0.1,0.1],[1,0]])
    with pytest.raises(TypeError):
        # define with invalid type for alphal 
        modpathrw.ModpathRWDsp(mp, alphal=None )
    with pytest.raises(TypeError):
        # define with invalid type for alphat 
        modpathrw.ModpathRWDsp(mp, alphat='alphat' )
    with pytest.raises(TypeError):
        # define with invalid type for dmeff 
        modpathrw.ModpathRWDsp(mp, dmeff='dmeff' )
    with pytest.raises(TypeError):
        # define with invalid a list with invalid interior type for dmeff 
        modpathrw.ModpathRWDsp(mp, dmeff=['dmeff'] )
    with pytest.raises(Exception):
        # define with inconsistent shape for alphalh array 
        modpathrw.ModpathRWDsp(mp, modelkind='axi', alphalh=[[0.1,0.1],[1,0]])
    with pytest.raises(Exception):
        # define with inconsistent shape for alphalv array 
        modpathrw.ModpathRWDsp(mp, modelkind='axi', alphalv=[[0.1,0.1],[1,0]])
    with pytest.raises(Exception):
        # define with inconsistent shape for alphath array 
        modpathrw.ModpathRWDsp(mp, modelkind='axi', alphath=[[0.1,0.1],[1,0]])
    with pytest.raises(Exception):
        # define with inconsistent shape for alphath array 
        modpathrw.ModpathRWDsp(mp, modelkind='axi', alphatv=[[0.1,0.1],[1,0]])
    with pytest.raises(TypeError):
        # define with invalid type for alphalh 
        modpathrw.ModpathRWDsp(mp, modelkind='axi', alphalh=None )
    with pytest.raises(TypeError):
        # define with invalid type for alphalhv
        modpathrw.ModpathRWDsp(mp, modelkind='axi', alphalv=None )
    with pytest.raises(TypeError):
        # define with invalid type for alphath
        modpathrw.ModpathRWDsp(mp, modelkind='axi', alphath=None )
    with pytest.raises(TypeError):
        # define with invalid type for alphatv 
        modpathrw.ModpathRWDsp(mp, modelkind='axi', alphavv=None )


    # Pass the dispersion parameters from the base case
    modpathrw.ModpathRWDsp(
        mp,
        alphal=MT3DP09Cases.alphal,
        alphat=MT3DP09Cases.alphat,
        dmeff=MT3DP09Cases.dmeff,
    )
    # Pass the dispersion parameters from the base case
    # for the axisymmetric form
    modpathrw.ModpathRWDsp(
        mp,
        modelkind = 'axi',
        alphalh   = MT3DP09Cases.alphal,
        alphalv   = MT3DP09Cases.alphal,
        alphath   = MT3DP09Cases.alphat,
        alphatv   = MT3DP09Cases.alphat,
        dmeff=MT3DP09Cases.dmeff,
    )

    dsp = modpathrw.ModpathRWDsp(mp)
    pkgs = mp.get_package_list() 
    assert dsp._ftype() in pkgs, (
            f"DSP package was not found in ModpathRW object"
        )
    # Pass the dispersion parameters from the base case
    modpathrw.ModpathRWDsp(
        mp,
        alphal=MT3DP09Cases.alphal,
        alphat=MT3DP09Cases.alphat,
        dmeff=MT3DP09Cases.dmeff,
    )

    with pytest.raises(Exception):
        # Try to write without sim, rwopts and bas package
        mp.write_input(check=True)

    # bas 
    modpathrw.ModpathRWBas(mp,porosity=MT3DP09Cases.porosity)

    # rwopts
    modpathrw.ModpathRWOpts(mp,dimensionsmask=[1,1,0])
    
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

    # Try to write ( checking consistency )
    mp.write_input(check=True)
