'''
Transport near low-permeability zone, backward.

Use a high-resolution flux concentration series of the extraction 
well as input for backward tracking of particles. Concentration data
is given as a timeseries aux variable in the flow model. 

References:
    - Zheng, C. & Wang, P., 1999, MT3DMS: a modular three-dimensional multispecies transport
      model for simulation of advection, dispersion, and chemical reactions of contaminants
      in groundwater systems; documentation and user’s guide.
    - https://modflow6-examples.readthedocs.io/en/master/_examples/ex-gwt-mt3dms-p09.html
'''

import os
import numpy as np 
import flopy
from flopy.utils.flopy_io import loadtxt
from flopyrw import modpathrw

# Model parameters #

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
k1 = 1.474e-4   # main permeability ($m/sec$)
k2 = 1.474e-7   # low permeability  ($m/sec$)
hk = k1 * np.ones((nlay, nrow, ncol), dtype=float)
hk[:, 5:8, 1:8] = k2 # Low permeability zone (delc=delr=100 [m])
icelltype = 0
idomain = np.ones((nlay, nrow, ncol), dtype=int)
# wells
qinjwell = 0.001
qextwell = -0.0189
cinjwell = 57.87
czero    = 0.0
injwell  = [0, 3, 6]
extwell  = [0, 10, 6]
# a single stress period with multiple time steps
nsteps  = 730
tlength = 2.0*365*86400
stress_periods = [ 
    {    
        'id': 0, 
        'length'       : tlength, 
        'n_time_steps' : nsteps, 
        'ts_multiplier': 1, 
        'steady_state' : False, 
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


# Configure the MODFLOW-6 model #

gwfname = "mf6model"
simname = "mprwexsim"
simws   = os.path.join(os.getcwd(),simname)

# flopy mf6 sim 
sim = flopy.mf6.MFSimulation(
    sim_name=simname, exe_name="mf6", version="mf6", sim_ws=simws
)

# tdis package
perioddata = []
for sp in stress_periods:
    perioddata.append( [sp['length'], sp['n_time_steps'], sp['ts_multiplier']] )
flopy.mf6.ModflowTdis(
    sim,
    nper=len(stress_periods),
    perioddata=perioddata,
    time_units=time_units, 
)

# gwf package
gwf = flopy.mf6.ModflowGwf(
    sim,
    modelname=gwfname,
    save_flows=True,
    model_nam_file="{}.nam".format(gwfname), 
)

# ims package
imsconfig = {
    'print_option'       : "SUMMARY",
    'outer_dvclose'      : hclose,
    'outer_maximum'      : nouter,
    'under_relaxation'   : "NONE",
    'inner_maximum'      : ninner,
    'inner_dvclose'      : hclose,
    'rcloserecord'       : rclose,
    'linear_acceleration': "CG",
    'scaling_method'     : "NONE",
    'reordering_method'  : "NONE",
    'relaxation_factor'  : relax,
    'filename'           : "{}.ims".format(gwfname),
}
imsgwf = flopy.mf6.ModflowIms(
        sim,
        **imsconfig
    )
sim.register_ims_package(imsgwf, [gwf.name]) 

# dis package
disconfig = { 
    'length_units': length_units, 
    'nlay'        : nlay,
    'nrow'        : nrow,
    'ncol'        : ncol,
    'delr'        : delr,
    'delc'        : delc,
    'top'         : top,
    'botm'        : botm,
    'idomain'     : idomain,
    'filename'    : "{}.dis".format(gwfname),
}
flopy.mf6.ModflowGwfdis(
    gwf,
    **disconfig
)

# ic package
flopy.mf6.ModflowGwfic(
    gwf,
    strt=chdh, 
    filename="{}.ic".format(gwfname)
)

# npf package
flopy.mf6.ModflowGwfnpf(  
    gwf,
    k=hk,
    k33=hk,
    save_flows=True,
    save_specific_discharge=True,
    filename="{}.npf".format(gwfname),
)

# chd package
xc = gwf.modelgrid.xcellcenters
chdspd = []
for j in range(ncol):
    # Loop through the top & bottom sides.
    #               l,  r, c,  head, conc
    chdspd.append([(0,  0, j), chdh, 0.0])  # Top boundary
    hd = 20.0 + (xc[-1, j] - xc[-1, 0]) * 2.5 / 100
    chdspd.append([(0, 17, j),    hd, 0.0])  # Bottom boundary

chdspd = {0: chdspd}
flopy.mf6.ModflowGwfchd(
    gwf,
    stress_period_data=chdspd,
    save_flows=True,
    maxbound=len(chdspd),
    auxiliary="CONCENTRATION",
    pname="CHD-1",
    filename="{}.chd".format(gwfname),
)

# wel package
# input cflux as timeseries.
# concentration is given with daily resolution (86400 s)
# until 2 years.
timeseries = {}
cfluxfname = 'p09mt3d_wout_cflux_highres.obs'
dtype = np.dtype(
    [
        ("tid"   , np.int32  ),
        ("time"  , np.float32),
        ("qsink" , np.float32),
        ("chist" , np.float32),
        ("cgpkde", np.float32),
    ]
)
cfluxdata = loadtxt( os.path.join( 'data', cfluxfname ), dtype=dtype, skiprows=0 )
times     = np.arange( 0.0, tlength+86400.0, 86400.0 )
cflux     = np.zeros( shape=times.shape )
cflux[ cfluxdata['tid']-1 ] = cfluxdata['cgpkde']
timeseries['timeseries'] = [ (t,cflux[it]) for it, t in enumerate( times ) ]
timeseries['time_series_namerecord']     = ['CFLUX']
timeseries['interpolation_methodrecord'] = ['LINEAR']

# first stress period
welsp1 = []
welsp1.append( # Injection well
    # (k, i, j), flow, conc
    [
        tuple(injwell), 
        qinjwell, 
        czero, 
    ]
) 
welsp1.append( # Extraction well
    [
        tuple(extwell), 
        qextwell, 
        'CFLUX', 
    ]
) 
welspd = {0: welsp1}
flopy.mf6.ModflowGwfwel(
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
flopy.mf6.ModflowGwfoc(
    gwf,
    head_filerecord="{}.hds".format(gwfname),
    budget_filerecord="{}.bud".format(gwfname),
    headprintrecord=[
        ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
    ],
    saverecord =[("HEAD", "ALL"), ("BUDGET", "ALL")],
    printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
)
# write and run 
sim.write_simulation()
sim.run_simulation()


# Configure the MODPATH-RW model #

# modpath-rw
mp = modpathrw.ModpathRW(
        modelname='mprwsim',
        flowmodel=gwf,
        model_ws =simws,
    )

# bas package
modpathrw.ModpathRWBas(mp,porosity=porosity)

# rwopts package
modpathrw.ModpathRWOpts(
    mp,
    timestep='min',
    ctdisp  =0.1,
    courant =0.1,
    dimensionsmask=[1,1,0],
)

# src package
sources = [
    (
        "WEL-1",
        [
            ["CONCENTRATION", 100.0, (8,8,1)],
        ],
    ),
]
src = modpathrw.ModpathRWSrc(
    mp,
    inputformat='aux', 
    sources=sources,
)

# dsp package
modpathrw.ModpathRWDsp(
    mp,
    alphal=alphal,
    alphat=alphat,
    dmeff =dmeff,
)

# obs package
# Note: the output time is modpath tracking time, 
# for interpretation in the context of the flow
# model it should be considered from the reference
# end time of modflow sim.
obs = modpathrw.ModpathRWObs(
    mp,
    kind='flux',
    cells=[tuple(injwell)],
    filename='obs_win_backward'
)

# sim
simconfig = {
    'simulationtype'         : 'rwendpoint', 
    'trackingdirection'      : 'backward',
    'weaksinkoption'         : 'pass_through',# particles pass at the extraction well (original, forward model) 
    'weaksourceoption'       : 'stop_at',     # particles stop at the injection well  (original, forward model)
    'referencetime'          : 2.0*365*86400, # modflow sim time is reversed
    'stoptimeoption'         : 'total',
    'endpointoutputoption'   : 2,
    'particlesmassoption'    : 0,
    'speciesdispersionoption': 0,
}
mprwsim = modpathrw.ModpathRWSim(
    mp, 
    **simconfig
)

# write and run 
mp.write_input()
mp.run_model(silent=False,report=True)
