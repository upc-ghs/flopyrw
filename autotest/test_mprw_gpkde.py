'''
Tests for the ModpathRWGpkde class
'''

import pytest
import numpy as np
from flopyrw import modpathrw
from autotest.test_mprw_mt3dp09_cases import MT3DP09Cases

def test_mprw_gpkde_input_mf6(function_tmpdir):
    '''
    Verifies the input for the GPKDE class 
    '''

    # get the mf6 case
    # brings WEL-1 and CHD-1 with aux CONCENTRATION
    gwf = MT3DP09Cases.mf6(function_tmpdir) 

    # modpath-rw
    mp = modpathrw.ModpathRW(
            modelname='mprwsim',
            flowmodel=gwf,
            model_ws =function_tmpdir,
        )
 
    # gpkde
    gpkde = modpathrw.ModpathRWGpkde(mp)

    # verify domainsize #
    #-------------------#
    xmin, xmax, ymin, ymax = gwf.modelgrid.extent
    zmax = np.max( gwf.modelgrid.top  )
    zmin = np.min( gwf.modelgrid.botm )

    # it should default to grid discretization params
    assert gpkde.domainsize[0] == abs(xmax-xmin), (
            f"gpkde: domainsize in x direction is not "
            f"consistent with expected grid dimensions"
        )
    assert gpkde.domainsize[1] == abs(ymax-ymin), (
            f"gpkde: domainsize in y direction is not "
            f"consistent with expected grid dimensions"
        )
    assert gpkde.domainsize[2] == abs(zmax-zmin), (
            f"gpkde: domainsize in z direction is not "
            f"consistent with grid discretization"
        )

    # verify binsize #
    #----------------#
    # it should default to grid discretization params
    assert MT3DP09Cases.delr == gpkde.binsize[0], (
            f"gpkde: binsize in x direction is not "
            f"consistent with grid discretization"
        )
    assert MT3DP09Cases.delc == gpkde.binsize[1], (
            f"gpkde: binsize in y direction is not "
            f"consistent with grid discretization"
        )
    assert MT3DP09Cases.delz == gpkde.binsize[2], (
            f"gpkde: binsize in z direction is not "
            f"consistent with grid discretization"
        )

    # verify domainorigin #
    #---------------------#
    assert gpkde.domainorigin[0] == gwf.modelgrid.xoffset , (
            f"gpkde: domainorigin in x direction is not "
            f"consistent with grid discretization"
        )
    assert gpkde.domainorigin[1] == gwf.modelgrid.yoffset , (
            f"gpkde: binsize in y direction is not "
            f"consistent with grid discretization"
        )
    assert gpkde.domainorigin[2] == np.min(gwf.modelgrid.botm), (
            f"gpkde: binsize in z direction is not "
            f"consistent with grid discretization"
        )


    # define consistent package 
    gpkde = modpathrw.ModpathRWGpkde(mp)

    # verify assignment to the main model
    pkgs = mp.get_package_list()
    assert gpkde._ftype() in pkgs, (
            f"GPKDE package was not found in ModpathRW object"
        )

    # and write (without checking model consistency, check=False by default)
    mp.write_input()



def test_mprw_gpkde_input_mf6disv(function_tmpdir):
    '''
    Verifies the input for the GPKDE class with disv grid 
    '''

    # get the mf6 case
    # brings WEL-1 and CHD-1 with aux CONCENTRATION
    gwf = MT3DP09Cases.mf6disv(function_tmpdir) 

    # modpath-rw
    mp = modpathrw.ModpathRW(
            modelname='mprwsim',
            flowmodel=gwf,
            model_ws =function_tmpdir,
        )
 
    # gpkde
    gpkde = modpathrw.ModpathRWGpkde(mp)

    # verify domainsize #
    #-------------------#
    xmin, xmax, ymin, ymax = gwf.modelgrid.extent
    zmax = np.max( gwf.modelgrid.top  )
    zmin = np.min( gwf.modelgrid.botm )

    # for a dis grid, without bin size, 
    # it should default to grid discretization params
    assert gpkde.domainsize[0] == abs(xmax-xmin), (
            f"gpkde: domainsize in x direction is not "
            f"consistent with expected grid dimensions"
        )
    assert gpkde.domainsize[1] == abs(ymax-ymin), (
            f"gpkde: domainsize in y direction is not "
            f"consistent with expected grid dimensions"
        )
    assert gpkde.domainsize[2] == abs(zmax-zmin), (
            f"gpkde: domainsize in z direction is not "
            f"consistent with grid discretization"
        )

    # verify binsize #
    #----------------#
    # it should default to grid discretization params
    assert MT3DP09Cases.delr == gpkde.binsize[0], (
            f"gpkde: binsize in x direction is not "
            f"consistent with grid discretization"
        )
    assert MT3DP09Cases.delc == gpkde.binsize[1], (
            f"gpkde: binsize in y direction is not "
            f"consistent with grid discretization"
        )
    assert MT3DP09Cases.delz == gpkde.binsize[2], (
            f"gpkde: binsize in z direction is not "
            f"consistent with grid discretization"
        )

    # verify domainorigin #
    #---------------------#
    assert gpkde.domainorigin[0] == gwf.modelgrid.xoffset , (
            f"gpkde: domainorigin in x direction is not "
            f"consistent with grid discretization"
        )
    assert gpkde.domainorigin[1] == gwf.modelgrid.yoffset , (
            f"gpkde: binsize in y direction is not "
            f"consistent with grid discretization"
        )
    assert gpkde.domainorigin[2] == np.min(gwf.modelgrid.botm), (
            f"gpkde: binsize in z direction is not "
            f"consistent with grid discretization"
        )

    # define consistent package 
    gpkde = modpathrw.ModpathRWGpkde(mp)

    # verify assignment to the main model
    pkgs = mp.get_package_list()
    assert gpkde._ftype() in pkgs, (
            f"GPKDE package was not found in ModpathRW object"
        )

    # and write (without checking model consistency, check=False by default)
    mp.write_input()


def test_mprw_gpkde_input_mf2005(function_tmpdir):
    '''
    Verifies the input for the GPKDE class 
    '''

    # get the mf2005 case
    # brings WEL with aux CONCENTRATION
    mf = MT3DP09Cases.mf2005(function_tmpdir)

    # modpath-rw
    mp = modpathrw.ModpathRW(
            modelname='mprwsimmf',
            flowmodel=mf,
            model_ws =function_tmpdir,
        )

    # gpkde
    gpkde = modpathrw.ModpathRWGpkde(mp)

    # verify domainsize #
    #-------------------#
    xmin, xmax, ymin, ymax = mf.modelgrid.extent
    zmax = np.max( mf.modelgrid.top  )
    zmin = np.min( mf.modelgrid.botm )

    # it should default to grid discretization params
    assert gpkde.domainsize[0] == abs(xmax-xmin), (
            f"gpkde: domainsize in x direction is not "
            f"consistent with expected grid dimensions"
        )
    assert gpkde.domainsize[1] == abs(ymax-ymin), (
            f"gpkde: domainsize in y direction is not "
            f"consistent with expected grid dimensions"
        )
    assert gpkde.domainsize[2] == abs(zmax-zmin), (
            f"gpkde: domainsize in z direction is not "
            f"consistent with grid discretization"
        )

    # verify binsize #
    #----------------#
    # it should default to grid discretization params
    assert MT3DP09Cases.delr == gpkde.binsize[0], (
            f"gpkde: binsize in x direction is not "
            f"consistent with grid discretization"
        )
    assert MT3DP09Cases.delc == gpkde.binsize[1], (
            f"gpkde: binsize in y direction is not "
            f"consistent with grid discretization"
        )
    assert MT3DP09Cases.delz == gpkde.binsize[2], (
            f"gpkde: binsize in z direction is not "
            f"consistent with grid discretization"
        )

    # verify domainorigin #
    #---------------------#
    assert gpkde.domainorigin[0] == mf.modelgrid.xoffset , (
            f"gpkde: domainorigin in x direction is not "
            f"consistent with grid discretization"
        )
    assert gpkde.domainorigin[1] == mf.modelgrid.yoffset , (
            f"gpkde: binsize in y direction is not "
            f"consistent with grid discretization"
        )
    assert gpkde.domainorigin[2] == np.min(mf.modelgrid.botm), (
            f"gpkde: binsize in z direction is not "
            f"consistent with grid discretization"
        )


    # define consistent package 
    gpkde = modpathrw.ModpathRWGpkde(mp)

    # verify assignment to the main model
    pkgs = mp.get_package_list()
    assert gpkde._ftype() in pkgs, (
            f"GPKDE package was not found in ModpathRW object"
        )

    # and write (without checking model consistency, check=False by default)
    mp.write_input()
