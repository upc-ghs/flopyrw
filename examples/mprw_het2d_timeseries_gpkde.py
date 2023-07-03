'''
Transport in two-dimensional heterogeneous aquifer.

Problem considers a rectangular initial condition of particles 
and transport through an heterogeneous system. 

Simulation considers smoothed reconstruction of concentrations for the 
final tracking time.
'''

import os
import flopy
import numpy as np
from flopyrw import modpathrw
from flopy.utils.flopy_io import loadtxt

# Model parameters #

hclose     = 1e-5
ninner     = 500 
nouter     = 100
nlay       = 1
nrow       = 50
ncol       = 250
hkvariance = 2.5
dtype = np.dtype(
    [
        ("lnk"   , np.float32 ),
    ]
)
lnkdata = loadtxt( os.path.join( 'data', 'lnk2ddata50x250.csv' ), dtype=dtype, skiprows=0 )
hkarray = np.exp( np.sqrt(hkvariance)*lnkdata['lnk'].reshape(nlay,nrow,ncol) )
stress_periods = [
        {
            'id'           : 0,
            'length'       : 100,
            'n_time_steps' : 1, 
            'ts_multiplier': 1,
            'steady_state' : True
        },
    ]
porosity = 0.35
dmeff    = 0
alphal   = 0.02
alphat   = 0.1*alphal


# Configure the MODFLOW-6 model #


simname   = 'mprwexsim'
modelname = 'mf6model'
simws     = os.path.join( os.getcwd(), simname )

# sim 
sim = flopy.mf6.MFSimulation(
        sim_name=simname,
        sim_ws  =simws,
        exe_name='mf6'
    )

# tdis package
perioddata = []
for sp in stress_periods:
    perioddata.append( [sp['length'], sp['n_time_steps'], sp['ts_multiplier']] )
tdis = flopy.mf6.ModflowTdis(
    sim,
    nper=len(stress_periods),
    perioddata=perioddata,
)

# ims package 
ims = flopy.mf6.ModflowIms(
        sim,
        pname               = 'ims',
        print_option        = 'SUMMARY',
        complexity          = 'MODERATE',
        inner_maximum       = ninner,
        inner_dvclose       = hclose, 
        linear_acceleration = 'BICGSTAB',
        scaling_method      = 'NONE',
        reordering_method   = 'NONE',
        relaxation_factor   = 0.97,
        outer_maximum       = nouter, 
        outer_dvclose       = hclose,
        no_ptcrecord        = ['ALL'], 
    )

# gwf package
gwf = flopy.mf6.ModflowGwf(
    sim,                                     
    modelname=modelname,
    save_flows=True,
    newtonoptions=['UNDER_RELAXATION'], # Newton
)

# dis package 
dis = flopy.mf6.ModflowGwfdis(
        gwf, 
        nlay=nlay, nrow=nrow, ncol=ncol, 
        delr=1, delc=1, top=1, botm=0, 
    )

# npf package
npf = flopy.mf6.ModflowGwfnpf(
    gwf,
    save_specific_discharge=True, 
    icelltype=0,
    k=hkarray,
)

# chd package
# unit mean gradient
hdin   = ncol+10
hdout  = 10
chdspd = []
for ir in range(nrow):
    chdspd.append( [ (0, ir,    0), hdin  ] )
    chdspd.append( [ (0, ir, ncol-1), hdout ] )
chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspd, 
        maxbound=len(chdspd),
    )

# ic package
# linear ic for each row
icarray= np.zeros(shape=(nlay,nrow,ncol), dtype=np.float64 )
xcoord = np.arange( 0, ncol, 1 ) 
for ir in range(nrow): 
    icarray[0,ir,:] = ( 1 - xcoord/ncol)*hdin + xcoord/ncol*hdout
ic = flopy.mf6.ModflowGwfic(
        gwf,
        strt=icarray,
    )

# oc package
budget_file = modelname + '.bud'
head_file   = modelname + '.hds'
flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=budget_file,
        head_filerecord  =head_file,
        saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')]
    )

# write and run 
sim.write_simulation()
success, mf6_output = sim.run_simulation(pause=False, report=True)
if not success:
    raise Exception('MF6 did not terminate normally !')


# Configure the MODPATH-RW model #


# modpath-rw
mp = modpathrw.ModpathRW(
        modelname = 'mprwsim',
        flowmodel = gwf,
        model_ws  = simws,
    )

# rwopts package
modpathrw.ModpathRWOpts(
        mp,
        timestep = 'min', 
        ctdisp   = 0.1,
        courant  = 0.1,
        dimensionsmask = [1,1,0],
    )

# basic package
modpathrw.ModpathRWBas(
        mp,
        porosity=porosity, 
    )

# ic package
xoinj = 10
dxinj = 10
yoinj = 10
dyinj = 30
concentration = np.zeros(shape=(nlay, nrow, ncol) ) 
concentration[0,yoinj:yoinj+dyinj,xoinj:xoinj+dxinj] = 1.0
modpathrw.ModpathRWIc(
        mp,
        concentration = concentration,
        particlesmass = 0.01, 
    )

# dsp package
modpathrw.ModpathRWDsp(
        mp,
        alphal=alphal, 
        alphat=alphat,
        dmeff = dmeff,
    )

# gpkde package 
gpkde = modpathrw.ModpathRWGpkde(
        mp, 
        skipinitialcondition = True, 
    )

# sim
simconfig = {
    'simulationtype'         : 'rwtimeseries', 
    'trackingdirection'      : 'forward',
    'referencetime'          : 0.0,
    'stoptimeoption'         : 'specified',
    'stoptime'               : 20.0,
    'timepointdata'          : [1,(20.0)],
    'timeseriesoutputoption' : 2, 
}
mprwsim = modpathrw.ModpathRWSim(
    mp, 
    **simconfig
)

# write and run 
mp.write_input()
mp.run_model(silent=False, report=True)

# get data to visualize the output
# Note.1: returns an array with the shape of the flowmodel grid. Up to this point 
#         works only for regular StructuredGrid, and reconstruction grid following
#         the flow-model discretization.
# Note.2: see flopyrw/modpathrw/mprwgpkde.py:ModpathRWGpkde:get_data for more details.
# Note.3: remember to plt.show() or savefig in order to visualize the figure.
cdata = gpkde.get_data()
pmv = flopy.plot.PlotMapView(gwf)
pmv.plot_array(cdata)

#import matplotlib.pyplot as plt
#plt.show()
