'''
'''
import numpy as np

from flopy import mf6

#from flopy.modflow import (
#    Modflow,
#    ModflowBas,
#    ModflowDis,
#    ModflowLpf,
#    ModflowOc,
#    ModflowPcg,
#    ModflowRch,
#    ModflowRiv,
#    ModflowWel,
#)
#from flopy.modpath import (
#    Modpath7,
#    Modpath7Bas,
#    Modpath7Sim,
#    ParticleData,
#    ParticleGroup,
#)


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


    ## TRANSPORT
    #porosity        = 0.3  # Porosity
    #sconc           = 0.0
    #al              = 20.0  # Longitudinal dispersivity ($m$)
    #trpt            = 0.2  # Ratio of horiz. transverse to longitudinal dispersivity ($m$)
    #ath1            = al * trpt
    #dmcoef          = 0.0  # m^2/s
    #sconc           = 0.0


    # Simple steady state stress periods
    stress_periods = [
            {
                'length'       : 1 * 365 * 86400, # 1 years, seconds
                'n_time_steps' : 1,
                'ts_multiplier': 1,
            },
            {
                'length'       : 1 * 365 * 86400, # 1 years, seconds
                'n_time_steps' : 1,
                'ts_multiplier': 1,
            },
        ]

    # solver
    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-6, 1.0


    @staticmethod
    def mf6(function_tmpdir):
        """
        Configures the flow model with MODFLOW 6
        """

        # model name
        ws = function_tmpdir / "mf6"
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
        for j in np.arange(MT3DP09Cases.ncol):
            # Loop through the top & bottom sides.
            #               l,  r, c,   head, conc
            chdspd.append([(0,  0, j), 250.0, 0.0])  # Top boundary
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
        # First stress period
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
        # Second stress period
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


        ## Write the datasets
        #sim.write_simulation()
        ## Run the simulation
        #success, buff = sim.run_simulation()
        #assert success, "mf6 model did not run"

        return gwf 




        ## create modpath files
        #mp = Modpath7(
        #    modelname=f"{nm}_mp", flowmodel=gwf, exe_name="mp7", model_ws=ws
        #)
        #defaultiface6 = {"RCH": 6, "EVT": 6}
        #mpbas = Modpath7Bas(mp, porosity=0.1, defaultiface=defaultiface6)
        #mpsim = Modpath7Sim(
        #    mp,
        #    simulationtype="combined",
        #    trackingdirection="forward",
        #    weaksinkoption="pass_through",
        #    weaksourceoption="pass_through",
        #    budgetoutputoption="summary",
        #    budgetcellnumbers=[1049, 1259],
        #    traceparticledata=[1, 1000],
        #    referencetime=[0, 0, 0.0],
        #    stoptimeoption="extend",
        #    timepointdata=[500, 1000.0],
        #    zonedataoption="on",
        #    zones=Mp7Cases.zones,
        #    particlegroups=Mp7Cases.particlegroups,
        #)

        ## write modpath datasets
        #mp.write_input()

        ## run modpath
        #success, buff = mp.run_model()
        #assert success, f"mp7 model ({mp.name}) did not run"

        #return mp

    #@staticmethod
    #def mf2005(function_tmpdir):
    #    """
    #    MODPATH 7 example 1 for MODFLOW-2005
    #    """

    #    ws = function_tmpdir / "mf2005"
    #    nm = "ex01_mf2005"
    #    iu_cbc = 130
    #    m = Modflow(nm, model_ws=ws, exe_name="mf2005")
    #    ModflowDis(
    #        m,
    #        nlay=Mp7Cases.nlay,
    #        nrow=Mp7Cases.nrow,
    #        ncol=Mp7Cases.ncol,
    #        nper=Mp7Cases.nper,
    #        itmuni=4,
    #        lenuni=1,
    #        perlen=Mp7Cases.perlen,
    #        nstp=Mp7Cases.nstp,
    #        tsmult=Mp7Cases.tsmult,
    #        steady=True,
    #        delr=Mp7Cases.delr,
    #        delc=Mp7Cases.delc,
    #        top=Mp7Cases.top,
    #        botm=Mp7Cases.botm,
    #    )
    #    ModflowLpf(
    #        m,
    #        ipakcb=iu_cbc,
    #        laytyp=Mp7Cases.laytyp,
    #        hk=Mp7Cases.kh,
    #        vka=Mp7Cases.kv,
    #    )
    #    ModflowBas(m, ibound=1, strt=Mp7Cases.top)
    #    # recharge
    #    ModflowRch(m, ipakcb=iu_cbc, rech=Mp7Cases.rch, nrchop=1)
    #    # wel
    #    wd = [i for i in Mp7Cases.wel_loc] + [Mp7Cases.wel_q]
    #    ModflowWel(m, ipakcb=iu_cbc, stress_period_data={0: wd})
    #    # river
    #    rd = []
    #    for i in range(Mp7Cases.nrow):
    #        rd.append(
    #            [
    #                0,
    #                i,
    #                Mp7Cases.ncol - 1,
    #                Mp7Cases.riv_h,
    #                Mp7Cases.riv_c,
    #                Mp7Cases.riv_z,
    #            ]
    #        )
    #    ModflowRiv(m, ipakcb=iu_cbc, stress_period_data={0: rd})
    #    # output control
    #    ModflowOc(
    #        m,
    #        stress_period_data={
    #            (0, 0): ["save head", "save budget", "print head"]
    #        },
    #    )
    #    ModflowPcg(m, hclose=1e-6, rclose=1e-3, iter1=100, mxiter=50)

    #    m.write_input()

    #    success, buff = m.run_model()
    #    assert success, "mf2005 model did not run"

    #    # create modpath files
    #    mp = Modpath7(
    #        modelname=f"{nm}_mp", flowmodel=m, exe_name="mp7", model_ws=ws
    #    )
    #    defaultiface = {"RECHARGE": 6, "ET": 6}
    #    mpbas = Modpath7Bas(mp, porosity=0.1, defaultiface=defaultiface)
    #    mpsim = Modpath7Sim(
    #        mp,
    #        simulationtype="combined",
    #        trackingdirection="forward",
    #        weaksinkoption="pass_through",
    #        weaksourceoption="pass_through",
    #        budgetoutputoption="summary",
    #        budgetcellnumbers=[1049, 1259],
    #        traceparticledata=[1, 1000],
    #        referencetime=[0, 0, 0.0],
    #        stoptimeoption="extend",
    #        timepointdata=[500, 1000.0],
    #        zonedataoption="on",
    #        zones=Mp7Cases.zones,
    #        particlegroups=Mp7Cases.particlegroups,
    #    )

    #    # write modpath datasets
    #    mp.write_input()

    #    return mp


    def case_mf6(self, function_tmpdir):
        return MT3DP09Cases.mf6(function_tmpdir)


    #def case_mf2005(self, function_tmpdir):
    #    return Mp7Cases.mf2005(function_tmpdir)
