'''
Tests for the ModpathRWSrc class
'''

import pytest
from flopyrw import modpathrw
from autotest.test_mprw_mt3dp09_cases import MT3DP09Cases

def test_mprw_src_aux_input_mf6(function_tmpdir):
    '''
    Verifies the inputformat aux with different
    source combinations, while using a mf6 model.
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
   
    # src

    # simple checks #
    #---------------#
    with pytest.raises(TypeError):
        # pass an int to inputformat
        modpathrw.ModpathRWSrc(
            mp,
            inputformat=0
        )
    with pytest.raises(ValueError):
        # pass a None to sources with aux format 
        modpathrw.ModpathRWSrc(
            mp,
            inputformat='aux',
            sources=None
        )
    with pytest.raises(ValueError):
        # pass an invalid input format 
        modpathrw.ModpathRWSrc(
            mp,
            inputformat='auxasd',
        )

    # pass invalid sources #
    #----------------------#
    with pytest.raises(ValueError):
        # pass sources with invalid pkg name
        sources = [
            (
                "VERYWEL-1", "CONCENTRATION", 100.0 , ( 4, 4, 1) 
            ),
        ]
        modpathrw.ModpathRWSrc(
            mp,
            inputformat='aux',
            sources=sources,
        )
    with pytest.raises(ValueError):
        # pass sources with invalid aux var name
        sources = [
            (
                "WEL-1", "THECONCENTRATION", 100.0 , ( 4, 4, 1) 
            ),
        ]
        modpathrw.ModpathRWSrc(
            mp,
            inputformat='aux',
            sources=sources,
        )
    with pytest.raises(ValueError):
        # pass sources with zero mass 
        sources = [
            (
                "WEL-1", "THECONCENTRATION", 0.0 , ( 4, 4, 1) 
            ),
        ]
        modpathrw.ModpathRWSrc(
            mp,
            inputformat='aux',
            sources=sources,
        )
    with pytest.raises(ValueError):
        # pass sources with invalid template
        sources = [
            (
                "WEL-1", "CONCENTRATION", 100.0 , ( 0, 4, 1) 
            ),
        ]
        modpathrw.ModpathRWSrc(
            mp,
            inputformat='aux',
            sources=sources,
        )
    with pytest.raises(ValueError):
        # pass sources with invalid uniform template 
        sources = [
            (
                "WEL-1", "CONCENTRATION", 100.0 , 0
            ),
        ]
        modpathrw.ModpathRWSrc(
            mp,
            inputformat='aux',
            sources=sources,
        )
    with pytest.raises(Exception):
        # pass repeated sources
        sources = [
            (
                "WEL-1", "CONCENTRATION", 100.0 , (4,4,1)
            ),
            (
                "WEL-1", "CONCENTRATION", 100.0 , (4,4,1)
            ),
        ]
        modpathrw.ModpathRWSrc(
            mp,
            inputformat='aux',
            sources=sources,
        )

    # Define a valid src
    sources = [
        (
            "WEL-1", "CONCENTRATION", 100.0 , (4,4,1)
        ),
    ]
    src = modpathrw.ModpathRWSrc(
        mp,
        inputformat='aux', 
        sources=sources,
    )
    pkgs = mp.get_package_list() 
    assert src._ftype() in pkgs, (
            f"SRC package was not found in ModpathRW object"
        )

    with pytest.raises(Exception):
        # Try to write without a sim package
        mp.write_input()

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

    # Try to write ( not checking consistency of dsp and rwopts, check=False )
    mp.write_input()

    
def test_mprw_src_aux_input_mf2005(function_tmpdir):
    '''
    Verifies the inputformat aux with different
    source combinations, while using a mf6 model.
    '''

    # get the mf2005 case
    # brings WEL with aux CONCENTRATION
    mf = MT3DP09Cases.mf2005(function_tmpdir)

    # modpath-rw
    mp = modpathrw.ModpathRW(
            modelname='mprwsim',
            flowmodel=mf,
            model_ws =function_tmpdir,
        )

    # src  
    with pytest.raises(Exception):
        # pass repeated sources
        sources = [
            (
                "WEL",
                [
                    ["CONCENTRATION", 100.0, (4,4,1)],
                    ["CONCENTRATION", 100.0, (4,4,1)]
                ],
            ),
        ]
        modpathrw.ModpathRWSrc(
            mp,
            inputformat='aux', 
            sources=sources,
        )

    with pytest.raises(ValueError):
        # pass invalid aux variable name
        sources = [
            (
                "WEL",
                [
                    ["THECONCENTRATION", 100.0, (4,4,1)],
                ],
            ),
        ]
        modpathrw.ModpathRWSrc(
            mp,
            inputformat='aux', 
            sources=sources,
        )
    with pytest.raises(ValueError):
        # pass invalid package name
        sources = [
            (
                "THEWEL",
                [
                    ["CONCENTRATION", 100.0, (4,4,1)],
                ],
            ),
        ]
        modpathrw.ModpathRWSrc(
            mp,
            inputformat='aux', 
            sources=sources,
        )


    # Define a valid source and write
    sources = [
        (
            "WEL",
            [
                ["CONCENTRATION", 100.0, (4,4,1)],
            ],
        ),
    ]
    src = modpathrw.ModpathRWSrc(
        mp,
        inputformat='aux', 
        sources=sources,
    )
    pkgs = mp.get_package_list() 
    assert src._ftype() in pkgs, (
            f"SRC package was not found in ModpathRW object"
        )

    with pytest.raises(Exception):
        # Try to write without a sim package
        mp.write_input()

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

    # Try to write ( not checking consistency of dsp and rwopts, check=False )
    mp.write_input()
