'''
Configures flow models with mf6 and mf2005 for P09 in MT3DMS
'''

import os
import numpy as np
from flopy import mf6
from flopy import modflow
from modflow_devtools.markers import requires_exe, requires_pkg


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

    # transport
    alphal   = 20 
    alphat   = 4
    dmeff    = 0.0
    porosity = 0.3


    @staticmethod
    def mf6(function_tmpdir, write=False, run=False, backward=False):
        """
        Configures the flow model with MODFLOW 6 and aux variables

        The backward enables non-zero concentration at extraction well
        and zero for injection.
        """

        # model name
        ws = function_tmpdir
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
        # modifies concentrations if backward tracking
        welsp1 = []
        if not backward:
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
        else:
            welsp1.append( # Injection well
                # (k, i, j), flow, conc
                [
                    tuple(MT3DP09Cases.injwell), 
                    MT3DP09Cases.qinjwell, 
                    MT3DP09Cases.czero
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
                    0.025*MT3DP09Cases.cinjwell # something not so large
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
    @requires_exe("gridgen")
    @requires_pkg("shapely", "shapefile")
    def mf6disv(function_tmpdir, write=False, run=False):
        """
        Configures the flow model with MODFLOW 6 and 
        a disv grid discretization
        """
        from flopy.utils import gridgen
        from shapely.geometry import Polygon
        from copy import deepcopy

        # model name
        ws = function_tmpdir
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

        # build aux dis grid
        auxgwf = deepcopy(gwf)
        delr = MT3DP09Cases.delr
        delc = MT3DP09Cases.delc
        disconfig = { 
            'nlay'        : MT3DP09Cases.nlay,
            'nrow'        : MT3DP09Cases.nrow,
            'ncol'        : MT3DP09Cases.ncol,
            'delr'        : MT3DP09Cases.delr,
            'delc'        : MT3DP09Cases.delc,
            'top'         : MT3DP09Cases.top,
            'botm'        : MT3DP09Cases.botm,
        }
        mf6.ModflowGwfdis(
                auxgwf,
                **disconfig
            )
        grid = gridgen.Gridgen(auxgwf.modelgrid, model_ws=ws)
        grid.build(verbose=False)
        gridprops = grid.get_gridprops_disv()

        # disv package
        mf6.ModflowGwfdisv(
            gwf,
            **gridprops,
            filename = "{}.disv".format(gwfname),
        )
        ncells = gwf.modelgrid.ncpl

        # ic package
        mf6.ModflowGwfic(
            gwf,
            strt=MT3DP09Cases.chdh, 
            filename="{}.ic".format(gwfname)
        )

        # npf package
        hk = MT3DP09Cases.k1 * np.ones(shape=(ncells,), dtype=float)
        pol= Polygon([ # low permeability zone
                ( 1*delr , 10*delc+1 ),
                ( 8*delr , 10*delc+1 ),
                ( 8*delr , 13*delc-1 ),
                ( 1*delr , 13*delc-1 ),
            ])
        gcells  = grid.intersect( [ pol ], 'polygon', 0 )
        nodes   = gcells['nodenumber']
        hk[ nodes ] = MT3DP09Cases.k2
        mf6.ModflowGwfnpf(  
            gwf,
            k=hk,
            k33=hk,
            save_flows=True,
            save_specific_discharge=True,
            filename="{}.npf".format(gwfname),
        )

        # chd package
        pol= Polygon([ # upper boundary
                ( 0*delr-1 , 17*delc+1 ),
                ( 14*delr+1, 17*delc+1 ),
                ( 14*delr+1, 18*delc-1 ),
                ( 0*delr-1 , 18*delc-1 ),
            ])
        gcells = grid.intersect( [ pol ], 'polygon', 0 )
        nodes  = gcells['nodenumber']
        xc     = gwf.modelgrid.xcellcenters

        chdspd = []
        for inode, node in enumerate(nodes):
            # [ ( lay, node ), head, caux ]
            chdspd.append([(0, node), MT3DP09Cases.chdh, 0.0])  # top boundary

        pol= Polygon([ # lower boundary
                ( 0*delr-1 , 0*delc+1 ),
                ( 14*delr+1, 0*delc+1 ),
                ( 14*delr+1, 0.5*delc-1 ),
                ( 0*delr-1 , 0.5*delc-1 ),
            ])
        gcells = grid.intersect( [ pol ], 'polygon', 0 )
        nodes  = gcells['nodenumber']
        for inod, nod  in enumerate( nodes ):
            hd = 20.0 + ( xc[nod] - xc[nodes[0]] ) * 2.5 / 100
            chdspd.append([(0, nod), hd, 0.0])  # bottom boundary
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
        polinj = Polygon([ # injection well 
                ( 600+1 , 1400+1 ),
                ( 700-1 , 1400+1 ),
                ( 700-1 , 1500-1 ),
                ( 600+1 , 1500-1 ),
            ])
        injcells  = grid.intersect( [ polinj ], 'polygon', 0 )
        injnodes  = injcells['nodenumber']
        nnodesinj = len(injnodes)

        polext = Polygon([ # extraction well
                ( 600+1 , 700+1 ),
                ( 700-1 , 700+1 ),
                ( 700-1 , 800-1 ),
                ( 600+1 , 800-1 ),
            ])
        extcells  = grid.intersect( [ polext ], 'polygon', 0 )
        extnodes  = extcells['nodenumber']
        nnodesext = len(extnodes)

        # first stress period
        welsp1 = []
        for ino, no in enumerate( injnodes ):
            welsp1.append( # injection well
                [
                    (0,no.item()),
                    MT3DP09Cases.qinjwell/nnodesinj, 
                    MT3DP09Cases.cinjwell
                ]
            ) 
        for ino, no in enumerate( extnodes ):
            welsp1.append( # extraction well
                [
                    (0,no.item()),
                    MT3DP09Cases.qextwell/nnodesext, 
                    MT3DP09Cases.czero
                ]
            ) 
        # second stress period
        welsp2 = []
        for ino, no in enumerate( injnodes ):
            welsp2.append( # injection well
                [
                    (0,no.item()),
                    MT3DP09Cases.qinjwell/nnodesinj, 
                    MT3DP09Cases.czero
                ]
            ) 
        for ino, no in enumerate( extnodes ):
            welsp2.append( # extraction well
                [
                    (0,no.item()),
                    MT3DP09Cases.qextwell/nnodesext, 
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
            assert success, "MT3DP09Cases: mf6disv model did not run."


        return gwf 


    @staticmethod
    @requires_exe("gridgen")
    @requires_pkg("shapely", "shapefile")
    def mf6disvmultilayer(function_tmpdir, write=False, run=False):
        """
        Configures the flow model with MODFLOW 6 and 
        a disv grid discretization
        """
        from flopy.utils import gridgen
        from shapely.geometry import Polygon
        from copy import deepcopy

        # model name
        ws = function_tmpdir
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

        # build aux dis grid
        auxgwf = deepcopy(gwf)
        delr = MT3DP09Cases.delr
        delc = MT3DP09Cases.delc
        nlay = 2*MT3DP09Cases.nlay
        delz = abs(MT3DP09Cases.top - MT3DP09Cases.botm[0])/nlay
        top  = MT3DP09Cases.top
        botm = np.arange(MT3DP09Cases.top-delz, MT3DP09Cases.botm[0]-delz, -delz).tolist()
        disconfig = { 
            'nlay'        : nlay,
            'nrow'        : MT3DP09Cases.nrow,
            'ncol'        : MT3DP09Cases.ncol,
            'delr'        : MT3DP09Cases.delr,
            'delc'        : MT3DP09Cases.delc,
            'top'         : top,
            'botm'        : botm,
        }
        mf6.ModflowGwfdis(
                auxgwf,
                **disconfig
            )
        grid = gridgen.Gridgen(auxgwf.modelgrid, model_ws=ws)
        grid.build(verbose=False)
        gridprops = grid.get_gridprops_disv()

        # disv package
        mf6.ModflowGwfdisv(
            gwf,
            **gridprops,
            filename = "{}.disv".format(gwfname),
        )
        ncells = gwf.modelgrid.ncpl

        # ic package
        mf6.ModflowGwfic(
            gwf,
            strt=MT3DP09Cases.chdh, 
            filename="{}.ic".format(gwfname)
        )

        # npf package
        hk = MT3DP09Cases.k1 * np.ones(shape=(nlay,ncells), dtype=float)
        pol= Polygon([ # low permeability zone
                ( 1*delr , 10*delc+1 ),
                ( 8*delr , 10*delc+1 ),
                ( 8*delr , 13*delc-1 ),
                ( 1*delr , 13*delc-1 ),
            ])
        gcells  = grid.intersect( [ pol ], 'polygon', 0 )
        nodes   = gcells['nodenumber']
        hk[ :, nodes ] = MT3DP09Cases.k2
        mf6.ModflowGwfnpf(  
            gwf,
            k=hk,
            k33=hk,
            save_flows=True,
            save_specific_discharge=True,
            filename="{}.npf".format(gwfname),
        )

        # chd package
        pol= Polygon([ # upper boundary
                ( 0*delr-1 , 17*delc+1 ),
                ( 14*delr+1, 17*delc+1 ),
                ( 14*delr+1, 18*delc-1 ),
                ( 0*delr-1 , 18*delc-1 ),
            ])
        gcells = grid.intersect( [ pol ], 'polygon', 0 )
        nodes  = gcells['nodenumber']
        xc     = gwf.modelgrid.xcellcenters

        chdspd = []
        for inode, node in enumerate(nodes):
            for ilay in range(nlay):
                # [ ( lay, node ), head, caux ]
                chdspd.append([(ilay, node), MT3DP09Cases.chdh, 0.0])  # top boundary

        pol= Polygon([ # lower boundary
                ( 0*delr-1 , 0*delc+1 ),
                ( 14*delr+1, 0*delc+1 ),
                ( 14*delr+1, 0.5*delc-1 ),
                ( 0*delr-1 , 0.5*delc-1 ),
            ])
        gcells = grid.intersect( [ pol ], 'polygon', 0 )
        nodes  = gcells['nodenumber']
        for inod, nod  in enumerate( nodes ):
            hd = 20.0 + ( xc[nod] - xc[nodes[0]] ) * 2.5 / 100
            for ilay in range(nlay):
                chdspd.append([(ilay, nod), hd, 0.0])  # bottom boundary
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
        polinj = Polygon([ # injection well 
                ( 600+1 , 1400+1 ),
                ( 700-1 , 1400+1 ),
                ( 700-1 , 1500-1 ),
                ( 600+1 , 1500-1 ),
            ])
        injcells  = grid.intersect( [ polinj ], 'polygon', 0 )
        injnodes  = injcells['nodenumber']
        nnodesinj = len(injnodes)

        polext = Polygon([ # extraction well
                ( 600+1 , 700+1 ),
                ( 700-1 , 700+1 ),
                ( 700-1 , 800-1 ),
                ( 600+1 , 800-1 ),
            ])
        extcells  = grid.intersect( [ polext ], 'polygon', 0 )
        extnodes  = extcells['nodenumber']
        nnodesext = len(extnodes)

        # arrays storing all well cells
        allinjnodes = np.zeros( shape=(nlay*nnodesinj,2), dtype=np.int32 )
        allextnodes = np.zeros( shape=(nlay*nnodesext,2), dtype=np.int32 )

        # first stress period
        welsp1 = []
        for ino, no in enumerate( injnodes ):
            for ilay in range(nlay):
                allinjnodes[ino*nlay+ilay,0] = ilay 
                allinjnodes[ino*nlay+ilay,1] = no
                welsp1.append( # injection well
                    [
                        (ilay,no.item()),
                        MT3DP09Cases.qinjwell/(nnodesinj*nlay), 
                        MT3DP09Cases.cinjwell
                    ]
                ) 
        for ino, no in enumerate( extnodes ):
            for ilay in range(nlay):
                allextnodes[ino*nlay+ilay,0] = ilay 
                allextnodes[ino*nlay+ilay,1] = no
                welsp1.append( # extraction well
                    [
                        (ilay,no.item()),
                        MT3DP09Cases.qextwell/(nnodesext*nlay),
                        MT3DP09Cases.czero
                    ]
                ) 
        # second stress period
        welsp2 = []
        for ino, no in enumerate( injnodes ):
            for ilay in range(nlay):
                welsp2.append( # injection well
                    [
                        (ilay,no.item()),
                        MT3DP09Cases.qinjwell/(nnodesinj*nlay), 
                        MT3DP09Cases.czero
                    ]
                ) 
        for ino, no in enumerate( extnodes ):
            for ilay in range(nlay):
                welsp2.append( # extraction well
                    [
                        (ilay,no.item()),
                        MT3DP09Cases.qextwell/(nnodesext*nlay), 
                        MT3DP09Cases.czero
                    ]
                ) 
        welspd = {0: welsp1, 1: welsp2}
        
        # wel
        mf6.ModflowGwfwel(
            gwf,
            stress_period_data=welspd,
            print_input = True,
            print_flows = True,
            save_flows  = True,
            auxiliary   = ["CONCENTRATION"],
            pname       = "WEL-1",
        )

        # pass allinjnodes and allextnodes to the class as nodenumbers 
        # for later use
        multiindex = tuple(allinjnodes.T)
        MT3DP09Cases.linearinjnodes = np.ravel_multi_index( multiindex , gwf.modelgrid.shape ) 
        multiindex = tuple(allextnodes.T)
        MT3DP09Cases.linearextnodes = np.ravel_multi_index( multiindex , gwf.modelgrid.shape )

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
            assert success, "MT3DP09Cases: mf6disvmultilayer model did not run."


        return gwf 



    @staticmethod
    @requires_exe("gridgen")
    @requires_pkg("shapely", "shapefile")
    def mf6disvusg(function_tmpdir, write=False, run=False):
        """
        Configures the flow model with MODFLOW 6 and 
        a disv grid discretization, with refinement, unstructured.
        """
        from flopy.utils import gridgen
        from shapely.geometry import Polygon
        from copy import deepcopy

        # model name
        ws = function_tmpdir
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
            'linear_acceleration': "BICGSTAB",
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

        # build aux dis grid
        auxgwf = deepcopy(gwf)
        delr = MT3DP09Cases.delr
        delc = MT3DP09Cases.delc
        disconfig = { 
            'nlay'        : MT3DP09Cases.nlay,
            'nrow'        : MT3DP09Cases.nrow,
            'ncol'        : MT3DP09Cases.ncol,
            'delr'        : MT3DP09Cases.delr,
            'delc'        : MT3DP09Cases.delc,
            'top'         : MT3DP09Cases.top,
            'botm'        : MT3DP09Cases.botm,
        }
        mf6.ModflowGwfdis(
                auxgwf,
                **disconfig
            )
        grid = gridgen.Gridgen(auxgwf.modelgrid, model_ws=ws)

        # refinement
        delr = MT3DP09Cases.delr
        delc = MT3DP09Cases.delc
        # injection area, level 1
        pol = Polygon([
                ( 3*delr   , 13*delc ),
                ( 3*delr   , 17*delc ),
                ( 10*delr  , 17*delc ),
                ( 10*delr  , 13*delc ),
            ])
        grid.add_refinement_features( [pol], 'polygon', 1, range(1) )
        # injection area, level 2
        pol = Polygon([
                ( 5*delr , 13*delc ),
                ( 5*delr , 16*delc ),
                ( 8*delr , 16*delc ),
                ( 8*delr , 13*delc ),
            ])
        grid.add_refinement_features( [pol], 'polygon', 2, range(1))
        # extraction area, level 1
        pol = Polygon([
                ( 3*delr   , 5*delc ),
                ( 3*delr   , 10*delc ),
                ( 10*delr  , 10*delc ),
                ( 10*delr  , 5*delc ),
            ])
        grid.add_refinement_features( [pol], 'polygon', 1, range(1) )
        # extraction area, level 2
        pol = Polygon([
                ( 5*delr , 6*delc ),
                ( 5*delr , 9*delc ),
                ( 8*delr , 9*delc ),
                ( 8*delr , 6*delc ),
            ])
        grid.add_refinement_features( [pol], 'polygon', 2, range(1))
        # low-permeability, level 1 
        pol = Polygon([
                ( 0*delr  , 9*delc+1  ),
                ( 10*delr , 9*delc+1  ),
                ( 10*delr , 14*delc-1 ),
                ( 0*delr  , 14*delc-1 ),
            ])
        grid.add_refinement_features( [pol], 'polygon', 1, range(1))
        # low-permeability, level 2 
        pol = Polygon([
                ( 5*delr , 9*delc+1 ),
                ( 9*delr , 9*delc+1 ),
                ( 9*delr , 14*delc-1 ),
                ( 5*delr , 14*delc-1 ),
            ])
        grid.add_refinement_features( [pol], 'polygon', 2, range(1))

        # build
        grid.build(verbose=False)
        gridprops = grid.get_gridprops_disv()

        # disv package
        mf6.ModflowGwfdisv(
            gwf,
            **gridprops,
            filename = "{}.disv".format(gwfname),
        )
        ncells = gwf.modelgrid.ncpl

        # ic package
        mf6.ModflowGwfic(
            gwf,
            strt=MT3DP09Cases.chdh, 
            filename="{}.ic".format(gwfname)
        )

        # npf package
        hk = MT3DP09Cases.k1 * np.ones(shape=(ncells,), dtype=float)
        pol= Polygon([ # low permeability zone
                ( 1*delr , 10*delc+1 ),
                ( 8*delr , 10*delc+1 ),
                ( 8*delr , 13*delc-1 ),
                ( 1*delr , 13*delc-1 ),
            ])
        gcells  = grid.intersect( [ pol ], 'polygon', 0 )
        nodes   = gcells['nodenumber']
        hk[ nodes ] = MT3DP09Cases.k2
        mf6.ModflowGwfnpf(  
            gwf,
            k=hk,
            k33=hk,
            xt3doptions = 'xt3d',
            save_flows=True,
            save_specific_discharge=True,
            filename="{}.npf".format(gwfname),
        )

        # chd package
        pol= Polygon([ # upper boundary
                ( 0*delr-1 , 17*delc+1 ),
                ( 14*delr+1, 17*delc+1 ),
                ( 14*delr+1, 18*delc-1 ),
                ( 0*delr-1 , 18*delc-1 ),
            ])
        gcells = grid.intersect( [ pol ], 'polygon', 0 )
        nodes  = gcells['nodenumber']
        xc     = gwf.modelgrid.xcellcenters

        chdspd = []
        for inode, node in enumerate(nodes):
            # [ ( lay, node ), head, caux ]
            chdspd.append([(0, node), MT3DP09Cases.chdh, 0.0])  # top boundary

        pol= Polygon([ # lower boundary
                ( 0*delr-1 , 0*delc+1 ),
                ( 14*delr+1, 0*delc+1 ),
                ( 14*delr+1, 0.5*delc-1 ),
                ( 0*delr-1 , 0.5*delc-1 ),
            ])
        gcells = grid.intersect( [ pol ], 'polygon', 0 )
        nodes  = gcells['nodenumber']
        for inod, nod  in enumerate( nodes ):
            hd = 20.0 + ( xc[nod] - xc[nodes[0]] ) * 2.5 / 100
            chdspd.append([(0, nod), hd, 0.0])  # bottom boundary
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
        polinj = Polygon([ # injection well 
                ( 600+1 , 1400+1 ),
                ( 700-1 , 1400+1 ),
                ( 700-1 , 1500-1 ),
                ( 600+1 , 1500-1 ),
            ])
        injcells  = grid.intersect( [ polinj ], 'polygon', 0 )
        injnodes  = injcells['nodenumber']
        nnodesinj = len(injnodes)

        polext = Polygon([ # extraction well
                ( 600+1 , 700+1 ),
                ( 700-1 , 700+1 ),
                ( 700-1 , 800-1 ),
                ( 600+1 , 800-1 ),
            ])
        extcells  = grid.intersect( [ polext ], 'polygon', 0 )
        extnodes  = extcells['nodenumber']
        nnodesext = len(extnodes)
        
        # pass injnodes and extnodes to the class 
        # for later use
        MT3DP09Cases.usginjnodes = injnodes
        MT3DP09Cases.usgextnodes = extnodes

        # first stress period
        welsp1 = []
        for ino, no in enumerate( injnodes ):
            welsp1.append( # injection well
                [
                    (0,no.item()),
                    MT3DP09Cases.qinjwell/nnodesinj, 
                    MT3DP09Cases.cinjwell
                ]
            ) 
        for ino, no in enumerate( extnodes ):
            welsp1.append( # extraction well
                [
                    (0,no.item()),
                    MT3DP09Cases.qextwell/nnodesext, 
                    MT3DP09Cases.czero
                ]
            ) 
        # second stress period
        welsp2 = []
        for ino, no in enumerate( injnodes ):
            welsp2.append( # injection well
                [
                    (0,no.item()),
                    MT3DP09Cases.qinjwell/nnodesinj, 
                    MT3DP09Cases.czero
                ]
            ) 
        for ino, no in enumerate( extnodes ):
            welsp2.append( # extraction well
                [
                    (0,no.item()),
                    MT3DP09Cases.qextwell/nnodesext, 
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
            assert success, "MT3DP09Cases: mf6disvusg model did not run."

        return gwf 



    @staticmethod
    def mf2005(function_tmpdir, write=False, run=False):
        """
        Configures the flow model with MODFLOW-2005
        """
        from copy import copy

        # model name
        ws        = function_tmpdir
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
            stress_period_data=chdspd,
            ipakcb=53
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
        # note: mp7 looks for budgetfilename
        # from lpf package
        modflow.ModflowLpf(
            mf,
            hk=MT3DP09Cases.hk,
            laytyp=MT3DP09Cases.laytyp,
            ss=0, sy=0,
            ipakcb=53,
        )

        options=[ 'AUX', 'CONCENTRATION' ]
        dtype = np.dtype(
                [
                    ("k", int),
                    ("i", int),
                    ("j", int),
                    ("flux", np.float32),
                    ("concentration", np.float32),
                ]
            )  
       
        # first stress period
        injwell = copy( MT3DP09Cases.injwell )
        injwell.extend( [MT3DP09Cases.qinjwell] ) # flux
        injwell.extend( [MT3DP09Cases.cinjwell] ) # concentration

        extwell = copy( MT3DP09Cases.extwell )
        extwell.extend( [MT3DP09Cases.qextwell] ) # flux 
        extwell.extend( [MT3DP09Cases.czero   ] ) # concentration
        welsp1 = [ injwell, extwell ]

        # second stress period
        injwell = copy( MT3DP09Cases.injwell )
        injwell.extend( [MT3DP09Cases.qinjwell] ) # flux
        injwell.extend( [MT3DP09Cases.czero   ] ) # concentration

        extwell = copy( MT3DP09Cases.extwell )
        extwell.extend( [MT3DP09Cases.qextwell] ) # flux 
        extwell.extend( [MT3DP09Cases.czero   ] ) # concentration
        welsp2 = [ injwell, extwell ]

        welspd = {0:welsp1, 1:welsp2}

        # wel
        wel = modflow.ModflowWel(
                mf,
                stress_period_data=welspd,
                dtype=dtype,
                options=options,
                ipakcb=53,
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



    @staticmethod
    def mf6multiaux(function_tmpdir, write=False, run=False):
        """
        Configures the flow model with MODFLOW 6 including 
        multiple aux variables:
         - CONC0: same values of the original test
         - CONC1: CONC0
         - CONC2: 2*CONC0
        """

        # model name
        ws = function_tmpdir
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
            # Loop through the top & bottom sides
            # lay, row, col, head, aux0, aux1, aux2
            chdspd.append([(0,  0, j), MT3DP09Cases.chdh, 0.0, 0.0, 0.0])  # Top boundary
            hd = 20.0 + (xc[-1, j] - xc[-1, 0]) * 2.5 / 100
            chdspd.append([(0, 17, j),    hd, 0.0, 0.0, 0.0])  # Bottom boundary
        chdspd = {0: chdspd}
        mf6.ModflowGwfchd(
            gwf,
            stress_period_data=chdspd,
            save_flows=True,
            maxbound=len(chdspd),
            auxiliary=["CONC0","CONC1","CONC2"],
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
                MT3DP09Cases.cinjwell,   # CONC0
                MT3DP09Cases.cinjwell,   # CONC1 = CONC0
                2*MT3DP09Cases.cinjwell, # CONC2 = 2*CONC0
            ]
        ) 
        welsp1.append( # Extraction well
            [
                tuple(MT3DP09Cases.extwell), 
                MT3DP09Cases.qextwell, 
                MT3DP09Cases.czero,
                MT3DP09Cases.czero,
                MT3DP09Cases.czero,
            ]
        ) 
        # second stress period
        welsp2 = []
        welsp2.append( # Injection well
            [
                tuple(MT3DP09Cases.injwell), 
                MT3DP09Cases.qinjwell, 
                MT3DP09Cases.czero,
                MT3DP09Cases.czero,
                MT3DP09Cases.czero
            ]
        ) 
        welsp2.append( # Extraction well
            [
                tuple(MT3DP09Cases.extwell), 
                MT3DP09Cases.qextwell, 
                MT3DP09Cases.czero,
                MT3DP09Cases.czero,
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
            auxiliary   = ["CONC0","CONC1","CONC2"],
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
    def mf6timeseriesaux(function_tmpdir, write=False, run=False, backward=False):
        """
        Configures the flow model with MODFLOW 6 including 
        a timeseries flow rate for the injection well, together 
        with a timeseries concentration. 
        multiple aux variables:
         - CONCENTRATION is the aux variable name. 
        """

        # model name
        ws = function_tmpdir
        nm = "p09mf6"
        gwfname = 'gwf-'+nm

        # flopy mf6 sim 
        sim = mf6.MFSimulation(
            sim_name=nm, exe_name="mf6", version="mf6", sim_ws=ws
        )

        # specific stress periods
        nsteps = 365
        stress_periods = [ 
            {    
                'id': 0, 
                'length'       : 365 * 86400, # 1 year, seconds 
                'n_time_steps' : nsteps, 
                'ts_multiplier': 1, 
                'steady_state' : False, 
            }, 
            {   
                'id': 1, 
                'length'       : 365 * 86400, # 1 year, seconds 
                'n_time_steps' : nsteps, 
                'ts_multiplier': 1, 
                'steady_state' : False, 
            } 
        ] 
        tlength = stress_periods[0]['length'] 
        deltat  = tlength/nsteps 
        times   = np.arange(0,tlength+deltat,deltat) 

        # tdis package
        perioddata = []
        for sp in stress_periods:
            perioddata.append( [sp['length'], sp['n_time_steps'], sp['ts_multiplier']] )
        mf6.ModflowTdis(
            sim,
            nper=len(stress_periods),
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
            # Loop through the top & bottom sides
            # lay, row, col, head, aux0, aux1, aux2
            chdspd.append([(0,  0, j), MT3DP09Cases.chdh, 0.0])  # Top boundary
            hd = 20.0 + (xc[-1, j] - xc[-1, 0]) * 2.5 / 100
            chdspd.append([(0, 17, j),    hd, 0.0])  # Bottom boundary
        chdspd = {0: chdspd}
        mf6.ModflowGwfchd(
            gwf,
            stress_period_data=chdspd,
            save_flows=True,
            maxbound=len(chdspd),
            auxiliary=["CONCENTRATION"],
            pname="CHD-1",
            filename="{}.chd".format(nm),
        )

        # Well input as timeseries
        timeseries = {}
        waveperiod = 10 * 86400  # n days in seconds
        omega      = np.pi/waveperiod
        phase      = 0.333*np.pi
        qinj       = MT3DP09Cases.qinjwell + 0.25*MT3DP09Cases.qinjwell*np.sin(omega*times) # qinjwell plus oscillation 
        cinj       = MT3DP09Cases.cinjwell*np.ones(shape=(len(times),))                     # constant c
        waveperiod = 30 * 86400  # n days in seconds
        omega      = np.pi/waveperiod
        qout       = MT3DP09Cases.qextwell + 0.05*MT3DP09Cases.qextwell*np.sin(omega*times) # qextwell plus oscillation
        timeseries['timeseries'] = [ (t,qinj[it],cinj[it], qout[it], 0.025*cinj[it]) for it, t in enumerate( times ) ]
        timeseries['time_series_namerecord']     = ['QIN', 'CIN', 'QOUT', 'CINOUT' ]
        timeseries['interpolation_methodrecord'] = ['STEPWISE','STEPWISE', 'STEPWISE', 'STEPWISE']
      
        # wel package
        # first stress period
        welsp1 = []
        if not backward:
            welsp1.append( # Injection well
                # (k, i, j), flow, conc
                [
                    tuple(MT3DP09Cases.injwell), 
                    'QIN', 
                    'CIN' 
                ]
            ) 
            welsp1.append( # Extraction well
                [
                    tuple(MT3DP09Cases.extwell), 
                    'QOUT', 
                    MT3DP09Cases.czero,
                ]
            ) 
            # second stress period
            welsp2 = []
            welsp2.append( # Injection well
                [
                    tuple(MT3DP09Cases.injwell), 
                    MT3DP09Cases.qinjwell,
                    MT3DP09Cases.czero,
                ]
            ) 
            welsp2.append( # Extraction well
                [
                    tuple(MT3DP09Cases.extwell), 
                    MT3DP09Cases.qextwell, 
                    MT3DP09Cases.czero,
                ]
            ) 
        else:
            welsp1.append( # Injection well
                # (k, i, j), flow, conc
                [
                    tuple(MT3DP09Cases.injwell), 
                    MT3DP09Cases.qinjwell,
                    MT3DP09Cases.czero,
                ]
            ) 
            welsp1.append( # Extraction well
                [
                    tuple(MT3DP09Cases.extwell), 
                    MT3DP09Cases.qextwell, 
                    MT3DP09Cases.czero,
                ]
            ) 
            # second stress period
            welsp2 = []
            welsp2.append( # Injection well
                [
                    tuple(MT3DP09Cases.injwell), 
                    MT3DP09Cases.qinjwell,
                    MT3DP09Cases.czero,
                ]
            ) 
            welsp2.append( # Extraction well
                [
                    tuple(MT3DP09Cases.extwell), 
                    MT3DP09Cases.qextwell, 
                    MT3DP09Cases.cinjwell,
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
            timeseries  = timeseries,
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




@requires_exe("mf6")
def test_mprw_mt3dp09_mf6(function_tmpdir):
    '''
    Write and run the mf6 model
    '''

    MT3DP09Cases.mf6(function_tmpdir,write=True,run=True)


@requires_exe("mf6", "gridgen")
def test_mprw_mt3dp09_mf6disv(function_tmpdir):
    '''
    Write and run the mf6 model with disv grid
    '''

    MT3DP09Cases.mf6disv(function_tmpdir,write=True,run=True)


@requires_exe("mf6", "gridgen")
def test_mprw_mt3dp09_mf6disvmultilayer(function_tmpdir):
    '''
    Write and run the mf6 model with disv grid
    '''

    MT3DP09Cases.mf6disvmultilayer(function_tmpdir,write=True,run=True)


@requires_exe("mf6", "gridgen")
def test_mprw_mt3dp09_mf6disvusg(function_tmpdir):
    '''
    Write and run the mf6 model with disvusg grid
    '''

    MT3DP09Cases.mf6disvusg(function_tmpdir,write=True,run=True)


@requires_exe("mf6")
def test_mprw_mt3dp09_mf6timeseriesaux(function_tmpdir):
    '''
    Write and run the mf6 model with disvusg grid
    '''

    MT3DP09Cases.mf6timeseriesaux(function_tmpdir,write=True,run=True)


@requires_exe("mf2005")
def test_mprw_mt3dp09_mf2005(function_tmpdir):
    '''
    Write and run the mf2005 model
    '''

    MT3DP09Cases.mf2005(function_tmpdir,write=True,run=True)

