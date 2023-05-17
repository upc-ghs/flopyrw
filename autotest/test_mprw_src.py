'''
Tests for the ModpathRWSrc class
'''

import pytest
from modflow_devtools.markers import requires_exe

from autotest.test_mprw_p09mt3d_cases import MT3DP09Cases
from flopy import modpathrw


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
