'''
Transport through a two-dimensional heterogeneous aquifer.

Problem considers a rectangular initial condition of particles 
and transport through an heterogeneous system. 

Simulation considers smoothed reconstruction of concentrations for several 
tracking times, employing a database of kernels to speed up the frequent 
concentration reconstruction process. The routine concludes with a short 
video showing the saved snapshots.

For further information regarding frequent spatial concentration reconstruction refer to:
  - Pérez-Illanes R., Fernàndez-Garcia, D., 2024, MODPATH-RW: A Random Walk Particle Tracking Code for Solute Transport in Heterogeneous Aquifers, Groundwater 62, no. 4: 617-634, doi:10.1111/gwat.13390
  - Pérez-Illanes R., Fernàndez-Garcia, D., 2024, A General Purpose Parallel Fortran Code for Grid Projected Concentration Reconstruction from Multidimensional Particle Distributions, Environmental Modelling & Software, doi: 10.1016/j.envsoft.2024.106008
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
try:
    import pandas as pd
    use_pandas = True
except ImportError:
    use_pandas = False
lnkdata = loadtxt(
    os.path.join( 'data', 'lnk2ddata50x250.csv' ),
    dtype=dtype,
    skiprows=0 ,
    use_pandas=use_pandas
)
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


#--------------------------------#
# Configure the MODPATH-RW model #
#--------------------------------#


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
# Note.1: by default the GPKDE package will define itself trying to follow 
#         the flowmodel extension, cell sizes.
# Note.2: kerneldatabase=True is recommended for models performing frequent 
#         concentration reconstruction.
gpkde = modpathrw.ModpathRWGpkde(
        mp,
        skipinitialcondition  = True,
        gridallocformat       = 1,
        convergence           = 0.01,
        initialsmoothingformat= 1,
        binsizefactor         = 1,
        kerneldatabase        = True,
        minhd                 = 0.25,
        deltahd               = 0.05,
        maxhd                 = 15,
    )

# sim
simconfig = {
    'simulationtype'         : 'rwtimeseries', 
    'trackingdirection'      : 'forward',
    'referencetime'          : 0.0,
    'stoptimeoption'         : 'specified',
    'stoptime'               : 30.0,
    'timepointdata'          : [20, 1.5],
    'timeseriesoutputoption' : 2, 
}
mprwsim = modpathrw.ModpathRWSim(
    mp, 
    **simconfig
)

# write and run 
mp.write_input()
mp.run_model(silent=False, report=True)

# plotting
import matplotlib.pyplot as plt
import matplotlib.colors as colors 

# get_alldata to visualize the output
# Note.1: returns an array with shape (ntimes,nlay,nrow,col). Up to this point 
#         is verified to work for regular StructuredGrid. For consistency, domain definition
#         in GPKDE should follow the flowmodel grid (i.e., cell size, domain extent).
# Note.2: see flopyrw/modpathrw/mprwgpkde.py:ModpathRWGpkde:get_alldata for more details.
# Note.3: remember to plt.show() or savefig in order to visualize the figures.
cdata = gpkde.get_alldata()

# create a colornorm based on the initial condition
norm = colors.Normalize(
        vmin=np.min(concentration[~np.isnan(concentration)]),
        vmax=np.max(concentration[~np.isnan(concentration)])
    )

# plot for all times. instead of pause and show one 
# could save each figure for subsequent postprocess/video.
pmv = flopy.plot.PlotMapView(gwf)
for it, time in enumerate(gpkde.times):
    if it > 0:
        im.set_array(cdata[it,0,:,:])
        im.set_norm( norm )
        plt.pause(0.3)
    else:
        im  = pmv.plot_array(cdata[it,0,:,:])
        im.set_norm( norm )
        pmv.ax.set_xlabel('x[m]')
        pmv.ax.set_ylabel('y[m]')
        plt.colorbar(im,ax=pmv.ax)
        plt.show(block=False)
