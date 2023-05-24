'''
Configures flow models with mf6 and mf2005 for P09 in MT3DMS
'''

import os
import numpy as np
from flopy import mf6
from flopy import modflow
from modflow_devtools.markers import requires_exe


class MT3DP09Cases:
    '''
    Setup flow model cases based on Problem 9 from MT3DMS

    https://modflow6-examples.readthedocs.io/en/master/_examples/ex-gwt-mt3dms-p09.html
    '''

    # units
    length_units = "meters"
    time_units   = "seconds"

    # grid
    delz  = 10.0
    delr  = 100.0
    delc  = 100.0
    sizey = 1800
    sizex = 1400
    nlay  = 1               
    nrow  = int(sizey/delr) 
    ncol  = int(sizex/delc) 
    top   = 0.0  
    botm  = [top - delz]
    chdh  = 250

    # flow
    k1      = 1.474e-4   # main permeability ($m/sec$)
    k2      = 1.474e-7   # low permeability  ($m/sec$)
    hk      = k1 * np.ones((nlay, nrow, ncol), dtype=float)
    hk[:, 5:8, 1:8] = k2 # Low permeability zone (delc=delr=100 [m])
    laytyp  = icelltype = 0
    idomain = np.ones((nlay, nrow, ncol), dtype=int)

    # wells
    qinjwell = 0.001
    qextwell = -0.0189
    cinjwell = 57.87
    czero    = 0.0
    injwell  = [0, 3, 6]
    extwell  = [0, 10, 6]

    # Simple steady state stress periods
    stress_periods = [
            {
                'length'       : 1 * 365 * 86400, # 1 years, seconds
                'n_time_steps' : 1,
                'ts_multiplier': 1,
                'steady'       : True,
            },
            {
                'length'       : 1 * 365 * 86400, # 1 years, seconds
                'n_time_steps' : 1,
                'ts_multiplier': 1,
                'steady'       : True,
            },
        ]

    # solver
    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-6, 1.0


    @staticmethod
    def mf6(function_tmpdir, write=False, run=False):
        """
        Configures the flow model with MODFLOW 6
        """

        # model name
        ws = os.path.join( function_tmpdir , "mf6")
        nm = "p09mf6"
        gwfname = 'gwf-'+nm

        # flopy mf6 sim 
        sim = mf6.MFSimulation(
            sim_name=nm, exe_name="mf6", version="mf6", sim_ws=ws
        )

        # tdis package
        perioddata = []
        for sp in MT3DP09Cases.stress_periods:
            perioddata.append( [sp['length'], sp['n_time_steps'], sp['ts_multiplier']] )
        mf6.ModflowTdis(
            sim,
            nper=len(MT3DP09Cases.stress_periods),
            perioddata=perioddata,
            time_units=MT3DP09Cases.time_units, 
        )

        # gwf package
        gwf = mf6.ModflowGwf(
            sim,
            modelname=gwfname,
            save_flows=True,
            model_nam_file="{}.nam".format(gwfname), 
        )

        # ims package
        imsconfig = {
            'print_option'       : "SUMMARY",
            'outer_dvclose'      : MT3DP09Cases.hclose,
            'outer_maximum'      : MT3DP09Cases.nouter,
            'under_relaxation'   : "NONE",
            'inner_maximum'      : MT3DP09Cases.ninner,
            'inner_dvclose'      : MT3DP09Cases.hclose,
            'rcloserecord'       : MT3DP09Cases.rclose,
            'linear_acceleration': "CG",
            'scaling_method'     : "NONE",
            'reordering_method'  : "NONE",
            'relaxation_factor'  : MT3DP09Cases.relax,
            'filename'           : "{}.ims".format(gwfname),
        }
        imsgwf = mf6.ModflowIms(
                sim,
                **imsconfig
            )
        sim.register_ims_package(imsgwf, [gwf.name]) 

        # dis package
        disconfig = { 
            'length_units': MT3DP09Cases.length_units, 
            'nlay'        : MT3DP09Cases.nlay,
            'nrow'        : MT3DP09Cases.nrow,
            'ncol'        : MT3DP09Cases.ncol,
            'delr'        : MT3DP09Cases.delr,
            'delc'        : MT3DP09Cases.delc,
            'top'         : MT3DP09Cases.top,
            'botm'        : MT3DP09Cases.botm,
            'idomain'     : MT3DP09Cases.idomain,
            'filename'    : "{}.dis".format(gwfname),
        }
        mf6.ModflowGwfdis(
            gwf,
            **disconfig
        )

        # ic package
        mf6.ModflowGwfic(
            gwf,
            strt=MT3DP09Cases.chdh, 
            filename="{}.ic".format(gwfname)
        )

        # npf package
        mf6.ModflowGwfnpf(  
            gwf,
            k=MT3DP09Cases.hk,
            k33=MT3DP09Cases.hk,
            save_flows=True,
            save_specific_discharge=True,
            filename="{}.npf".format(gwfname),
        )

        # chd package
        xc = gwf.modelgrid.xcellcenters
        chdspd = []
        for j in range(MT3DP09Cases.ncol):
            # Loop through the top & bottom sides.
            #               l,  r, c,               head, conc
            chdspd.append([(0,  0, j), MT3DP09Cases.chdh, 0.0])  # Top boundary
            hd = 20.0 + (xc[-1, j] - xc[-1, 0]) * 2.5 / 100
            chdspd.append([(0, 17, j),    hd, 0.0])  # Bottom boundary
        chdspd = {0: chdspd}
        mf6.ModflowGwfchd(
            gwf,
            stress_period_data=chdspd,
            save_flows=True,
            maxbound=len(chdspd),
            auxiliary="CONCENTRATION",
            pname="CHD-1",
            filename="{}.chd".format(nm),
        )

        # wel package
        # first stress period
        welsp1 = []
        welsp1.append( # Injection well
            # (k, i, j), flow, conc
            [
                tuple(MT3DP09Cases.injwell), 
                MT3DP09Cases.qinjwell, 
                MT3DP09Cases.cinjwell
            ]
        ) 
        welsp1.append( # Extraction well
            [
                tuple(MT3DP09Cases.extwell), 
                MT3DP09Cases.qextwell, 
                MT3DP09Cases.czero
            ]
        ) 
        # second stress period
        welsp2 = []
        welsp2.append( # Injection well
            [
                tuple(MT3DP09Cases.injwell), 
                MT3DP09Cases.qinjwell, 
                MT3DP09Cases.czero
            ]
        ) 
        welsp2.append( # Extraction well
            [
                tuple(MT3DP09Cases.extwell), 
                MT3DP09Cases.qextwell, 
                MT3DP09Cases.czero
            ]
        ) 
        welspd = {0: welsp1, 1: welsp2}
        mf6.ModflowGwfwel(
            gwf,
            stress_period_data=welspd,
            print_input = True,
            print_flows = True,
            save_flows  = True,
            auxiliary   = ["CONCENTRATION"],
            pname       = "WEL-1",
        )

        # oc package
        mf6.ModflowGwfoc(
            gwf,
            head_filerecord="{}.hds".format(gwfname),
            budget_filerecord="{}.bud".format(gwfname),
            headprintrecord=[
                ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
            ],
            saverecord =[("HEAD", "ALL"), ("BUDGET", "ALL")],
            printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )

        if write:
            sim.write_simulation(silent=True)
        if run:
            success, buff = sim.run_simulation(silent=True)
            assert success, "MT3DP09Cases: mf6 model did not run."

        return gwf 


    @staticmethod
    def mf2005(function_tmpdir, write=False, run=False):
        """
        Configures the flow model with MODFLOW-2005
        """
        from copy import copy

        # model name
        ws        = os.path.join( function_tmpdir , "mf2005" )
        nm        = "p09mf"
        modelname = 'mf-'+nm

        # modflow
        mf = modflow.Modflow(
            modelname=modelname, model_ws=ws, exe_name="mf2005"
        )

        # Time variables
        nper   = len(MT3DP09Cases.stress_periods)
        perlen = []
        nstp   = []
        tsmult = []
        steady = []
        for sp in MT3DP09Cases.stress_periods:
            perlen.append( sp['length'] )
            nstp.append( sp['n_time_steps'] )
            tsmult.append( sp['ts_multiplier'] )
            steady.append( sp['steady'] )

        # dis package
        disconfig = { 
            'nlay'        : MT3DP09Cases.nlay,
            'nrow'        : MT3DP09Cases.nrow,
            'ncol'        : MT3DP09Cases.ncol,
            'delr'        : MT3DP09Cases.delr,
            'delc'        : MT3DP09Cases.delc,
            'top'         : MT3DP09Cases.top,
            'botm'        : MT3DP09Cases.botm,
            'itmuni'      : 1, # seconds
            'lenuni'      : 2, # meters
            'nper'        : nper,
            'perlen'      : perlen,
            'steady'      : steady,
            'nstp'        : nstp,
        }
        modflow.ModflowDis(
                mf,
                **disconfig 
            )

        # chd 
        xc = mf.modelgrid.xcellcenters
        chdspd = []
        for j in range(MT3DP09Cases.ncol):
            # Loop through the top & bottom sides.
            topcell = [0,0,j] # top boundary
            topcell.extend( 2*[MT3DP09Cases.chdh] )
            chdspd.append( topcell ) 
            hd = 20.0 + (xc[-1, j] - xc[-1, 0]) * 2.5 / 100
            botcell = [0,17,j]
            botcell.extend( 2*[hd] )
            chdspd.append(botcell) # Bottom boundary
        chdspd = {0: chdspd}
        modflow.ModflowChd(
            mf,
            stress_period_data=chdspd
        )

        # bas 
        ibound = np.ones(
                (
                    MT3DP09Cases.nlay,
                    MT3DP09Cases.nrow,
                    MT3DP09Cases.ncol
                ),
                dtype=np.int32
            )
        ibound[0, 0, :]  = -1
        ibound[0, -1, :] = -1
        strt = np.ones(
                (
                    MT3DP09Cases.nlay,
                    MT3DP09Cases.nrow,
                    MT3DP09Cases.ncol
                ),
                dtype=np.float32
            )*MT3DP09Cases.chdh
        for j in range(MT3DP09Cases.ncol):
            strt[0, -1, j] = 20.0 + (xc[-1, j] - xc[-1, 0]) * 2.5 / 100
        modflow.ModflowBas(
            mf,
            ibound=ibound,
            strt=strt,
        )

        # lpf
        modflow.ModflowLpf(
            mf,
            hk=MT3DP09Cases.hk,
            laytyp=MT3DP09Cases.laytyp,
            ss=0, sy=0,
        )

        # wel
        injwell = copy( MT3DP09Cases.injwell )
        injwell.extend( [MT3DP09Cases.qinjwell] )
        extwell = copy( MT3DP09Cases.extwell )
        extwell.extend( [MT3DP09Cases.qextwell] )
        welspd = {
            0: [ injwell, extwell ]
        }
        modflow.ModflowWel(
            mf,
            stress_period_data=welspd
        )

        # pcg
        modflow.ModflowPcg(mf)

        # oc
        modflow.ModflowOc(
            mf,
            stress_period_data={
                (0,0):['save head', 'save budget'],
                (1,0):['save head', 'save budget']
            }
        )
        
        if write:
            mf.write_input()
        if run:
            success, buff = mf.run_model(silent=True)
            assert success, "MT3DP09Cases: mf2005 model did not run."

        return mf


    def case_mf6(self, function_tmpdir, write=False, run=False):
        return MT3DP09Cases.mf6(function_tmpdir, write=write, run=run)

    def case_mf2005(self, function_tmpdir, write=False, run=False):
        return MT3DP09Cases.mf2005(function_tmpdir, write=write, run=run)


@requires_exe("mf6")
def test_mprw_mt3dp09_mf6(function_tmpdir):
    '''
    Write and run the mf6 model
    '''

    MT3DP09Cases.mf6(function_tmpdir,write=True,run=True)


@requires_exe("mf2005")
def test_mprw_mt3dp09_mf2005(function_tmpdir):
    '''
    Write and run the mf2005 model
    '''

    MT3DP09Cases.mf2005(function_tmpdir,write=True,run=True)
