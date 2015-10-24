from woo.core import *
from woo.dem import *
from minieigen import *
import woo, woo.models
import woo.pre.triax
import woo.pre.cylTriax
woo.master.usesApi=10102

S=woo.master.scene=Scene(fields=[DemField(gravity=(0,0,-9.8))])
#mat=woo.dem.FrictMat(young=1e7,ktDivKn=.2,tanPhi=.4);
mat=woo.dem.FrictMat(density=2150, id=-1, young=50e6, tanPhi=0.5, ktDivKn=0.25);

partMask=0b0001;
wallMask=0b0011;		#DemField.defaultBoundaryMask
loneMask=0b0010;		#DemField.defaultLoneMask 
S.dem.loneMask=loneMask;

model=woo.models.ContactModelSelector(
    name='linear',
    damping=0.4,
    mats=[mat]
    )

S.lab.wallMat=model.mats[0].deepcopy()
S.lab.wallMat.tanPhi=0

iniSize=(0.4,0.2,0.8);
box=AlignedBox3((0,0,0),iniSize);
#box=AlignedBox3((0,0,0),(2,2,2))
#gen=PsdSphereGenerator(psdPts=[(.1,0),(.1,.2),(.2,.9),(.3,1)])
gen=MinMaxSphereGenerator(dRange=(0.02,0.02));   
#ctrlStep=100
ctrlStep=100;

factoryKw=dict(
	generator=gen,
	materials=[mat],
	massRate=0,
	maxMass=-1,
	maxNum=-1,
	maxAttempts=5000,
        atMaxAttempts='dead',
	shooter=None,
	mask=partMask,
        label='inlet',
	collideExisting=False
	);

S.trackEnergy=True;
S.dem.par.add(Wall.makeBox(box=box,which=(1,1,1,1,1,0),mat=S.lab.wallMat));
#engine 1: create particles
S.engines=[
    BoxInlet(
        stepPeriod=ctrlStep,
        box=box,
	**factoryKw
        ),
#    PyRunner(100,'print(S.step,S.time,S.dt)'),
#    PyRunner(command='addPlotData(S)',stepPeriod=100),
#    PyRunner(100,'addPlotData(S)'),
#    PyRunner(command='checkUnbalanced(S)',realPeriod=2)
    ];

#S.plot.plots={'i':('unbalanced',None,'**S.energy')};
#S.plot.plot()
S.saveTmp();
#S.run();
#S.run(21000,True);
S.one();

print 'Number of nodes:',len(S.dem.nodes)
S.periodic=True
margin=0
S.cell.setBox((1+2*margin)*iniSize)
#S.periodic=True
#S.cell.setBox((1+2*margin)*pre.iniSize)
#S.dtSafety=pre.dtSafety


# parameters to be used in the engine 2
sigIso=-200e3;
maxRates=(2e-1,2e-1,1.0);
rateStep=0.01;
stopStrain=-0.3;
maxUnbalanced=0.1;
massFactor=0.2;
S.lab.relVol=1.0;
S.dtSafety=0.7;
#engine 2: compaction and shearing
S.engines=[
    woo.dem.WeirdTriaxControl(
	  goal=(sigIso,sigIso,sigIso),
	  maxStrainRate=(maxRates[0],maxRates[0],maxRates[0]),
	  relVol=S.lab.relVol,
	  stressMask=0b01111,
	  globUpdate=1,
	  maxUnbalanced=maxUnbalanced,
	  mass=massFactor*sum([n.dem.mass for n in S.dem.nodes]),
	  doneHook='compactionDone(S)',
      label='triax',
      absStressTol=1e4,
      relStressTol=1e-2
    ),
#    woo.core.PyRunner(20,'addPlotData_checkProgress(S)')
     woo.core.PyRunner(20,'import woo.pre.cylTriax;addPlotData_checkProgress(S)')
]+woo.utils.defaultEngines(model=model,dynDtPeriod=100)
        
S.lab.stage='compact'
S.lab._setWritable('stage')

def addPlotData_checkProgress(S):
    assert S.lab.stage in ('compact','triax')
#    import woo
    t=S.lab.triax

    sxx,syy,szz=t.stress.diagonal() 
    dotE=S.cell.gradV.diagonal()
    dotEMax=t.maxStrainRate
    # net volume
    vol=S.cell.volume*S.lab.relVol
    # current radial stress
    srr=.5*(sxx+syy) 
    # mean stress
    p=t.stress.diagonal().sum()/3.
    # deviatoric stress
    q=szz-.5*(sxx+syy) 
    qDivP=(q/p if p!=0 else float('nan'))
    SzzDivSxx=(szz/sxx if sxx!=0 else float('nan'))

    if S.lab.stage=='compact':
        ## t.strain is log(l/l0) for all components
        exx,eyy,ezz=t.strain 
        err=.5*(exx+eyy)
        # volumetric strain is not defined directly, and it is not needed either        
        eVol=float('nan')
    else:
        # triaxial phase:
        # only axial strain (ezz) and volumetric strain (eVol) are known
        #
        # set the initial volume, if not yet done
        if not hasattr(S.lab,'netVol0'): S.lab.netVol0=S.cell.volume*S.lab.relVol
        # axial strain is known; xy components irrelevant (inactive)
        ezz=t.strain[2] 
        # current volume / initial volume
        eVol=math.log(vol/S.lab.netVol0) 
        # radial strain
        err=.5*(eVol-ezz) 
        # undefined
        exx=eyy=float('nan') 

    # deviatoric strain in z direction
    eDev=ezz-(1/3.)*(2*err+ezz)

    S.plot.addData(
        unbalanced=woo.utils.unbalancedForce(),
        i=S.step,
        time=S.time,
        sxx=sxx,syy=syy,srr=.5*(sxx+syy),szz=szz,
        exx=exx,eyy=eyy,err=err,ezz=ezz,
        dotE=dotE,dotErr=.5*(dotE[0]+dotE[1]),
        dotEMax=dotEMax,
        dotEMax_z_neg=-dotEMax[2],
        eDev=eDev,eVol=eVol,
        vol=vol,
        p=p,q=q,qDivP=qDivP,
        SzzDivSxx=SzzDivSxx,
        
        isTriax=(1 if S.lab.stage=='triax' else 0), # to be able to filter data
        # parTanPhi=S.lab.partMat.tanPhi,
        #memTanPhi=S.lab.memMat.tanPhi,
        #suppTanPhi=S.lab.suppMat.tanPhi
        # save all available energy data
        #Etot=O.energy.total()#,**O.energy
    )

    if not S.plot.plots:
        S.plot.plots={
            'i':('unbalanced',None,'vol'),
            ' i':('sxx','syy','szz'),
            'i ':('exx','eyy','ezz','eVol'),
#            ' i ':('dotE_z','dotEMax_z'),
            'ezz':(('qDivP','g-'),None,('eVol','r-')),
            'time':('ezz',),
            ' ezz':('SzzDivSxx',),
            # energy plot
            #' i ':(O.energy.keys,None,'Etot'),
        };
        S.plot.xylabels={
            'i':('step','Stress(Pa)','Volume'),
            ' i':('step','Stress(Pa)',),
            'i ':('step', 'strains'),
 #   ' i ':('step','Strains','Strains'),
            'ezz':('Axial strain','Ratio','Volumetic strain'),
            'time':('time','Axial strain'),
        };
        S.plot.labels={
            'sxx':r'$\sigma_{xx}$',
            'syy':r'$\sigma_{yy}$',
            'szz':r'$\sigma_{zz}$',
            'srr':r'$\sigma_{rr}$',
            'surfLoad':r'$\sigma_{\rm hydro}$',
            'exx':r'$\varepsilon_{xx}$',
            'eyy':r'$\varepsilon_{yy}$',
            'ezz':r'$\varepsilon_{zz}$',
            'err':r'$\varepsilon_{rr}$',
            'eVol':r'$\varepsilon_{v}$',
            'vol':'volume',
            'eDev':r'$\varepsilon_d$',
            'qDivP':'$q/p$',
            'p':'$p$',
            'q':'$q$',
            'dotE_x':r'$\dot\varepsilon_{xx}$',
            'dotE_y':r'$\dot\varepsilon_{yy}$',
            'dotE_z':r'$\dot\varepsilon_{zz}$',
            'dotE_rr':r'$\dot\varepsilon_{rr}$',
            'dotEMax_z':r'$\dot\varepsilon_{zz}^{\rm max}$',
            'dotEMax_z_neg':r'$-\dot\varepsilon_{zz}^{\rm max}$'
        };

    ## adjust rate in the triaxial stage
    if S.lab.stage=='triax':
        t.maxStrainRate[2]=min(t.maxStrainRate[2]+rateStep*maxRates[1],maxRates[1])


def compactionDone(S):
    # if S.lab.compactMemoize:
    print 'Compaction done at step',S.step
    import woo
    t=S.lab.triax
    # set the current cell configuration to be the reference one
    S.cell.trsf=Matrix3.Identity
    S.cell.refHSize=S.cell.hSize
    S.cell.nextGradV=Matrix3.Zero # avoid spurious strain from the last gradV value
    ##  S.lab.leapfrog.damping=.7 # increase damping to a larger value
    t.stressMask=0b0001 # z is strain-controlled, y has zero strain, x is stress-controlled
    t.goal=(sigIso,0,stopStrain)
    # start with small rate, increases later
    t.maxStrainRate=(maxRates[2],maxRates[2],rateStep*maxRates[1]) 
    t.maxUnbalanced=10 # don't care about that
    t.doneHook='triaxDone(S)'

    # recover friction angle
    S.lab.partMat.tanPhi=S.model.mats[0].tanPhi
    # force update of contact parameters
    S.lab.contactLoop.updatePhys='once'


    try:
        import woo.gl
        woo.gl.Gl1_DemField.updateRefPos=True
    except ImportError: pass

    S.lab.stage='triax'

    if S.pre.saveFmt:
        out=S.pre.saveFmt.format(stage='compact',S=S,**(dict(S.tags)))
        print 'Saving to',out
        S.save(out)


def triaxDone(S):
    print 'Triaxial done at step',S.step
    if S.pre.saveFmt:
        out=S.pre.saveFmt.format(stage='done',S=S,**(dict(S.tags)))
        print 'Saving to',out
        S.save(out)
    S.stop()
    import woo.utils
    (repName,figs)=woo.utils.htmlReport(S,S.pre.reportFmt,'Triaxial test with rigid boundary',afterHead='',figures=[(None,f) for f in S.plot.plot(noShow=True,subPlots=False)],svgEmbed=True,show=True)
    woo.batch.writeResults(S,defaultDb='triax.hdf5',series=S.plot.data,postHooks=[woo.pre.cylTriax.plotBatchResults],simulationName='triax',report='file://'+repName)


