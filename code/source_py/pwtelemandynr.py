import os
from os.path import join, isfile, islink, getmtime
from cmath import exp
import array

import numpy as np

from ase.data import chemical_symbols
from ase.units import Bohr,Rydberg, fs,Hartree
from ase.io.cube import read_cube_data

class pwtelemandynr:
    def __init__(self, label='pwteleman',xc='LDA',
charge=None,
basis=None,
kxbox=36,
kybox=36,
kzbox=36,
kstate=100,
nspdw=50,
dx=-0.8,
dy=-0.8,
dz=-0.8,
iforce=0,
ipseudo=1,          
epswf=0.2,
e0dmp=2.0,
epsorc=1e-8 ,      
b2occ=0.4,
gamocc=10.0,
deocc=0.0,
osfac=1.0 ,
temp=0.0,
occmix=0.5,
epsoro=1e-6,    
isurf=0,
nk=100,
nc=100,
endcon=1e-5,
radjel=4e0,
surjel=1e0,    
itback=200,                            
bbeta=0.0,
gamma=0.0,
beta4=0.0,       
dpolx=0.0,
dpoly=0.0,
dpolz=0.0,      
iexcit=0,
irotat=0,
ispidi=0,         
phirot=0.0,
scaleclust=1.0,
scaleclustx=1.0,
scaleclusty=1.0,
scaleclustz=1.0,
iemomsrel=1,
shiftclustx=0.0,
shiftclusty=0.0,
shiftclustz=0.0,
rotclustx=0.0,
rotclusty=0.0,
rotclustz=0.0,
imob=1,
shiftwfx=0.0,
shiftwfy=0.0,
shiftwfz=0.0,
ispinsep=0,
phangle=0.0,
phphase=0.0,       
idebug=0,
ifreezekspot=0,
itof=0,
jescmask=0,
ishiftcmtoorigin=0,
jescmaskorb=0,
ishutdown=0,
jplotdensitydiff=0,
jplotdensity2d=0,
jplotdensitydiff2d=0,
nmptheta=2,
nmpphi=1,
jmp=0,
jovlp=100000000,
jnorms=0,
iscatterelectron=0,
jcharges=0,
iaddcluster=0,
iswforce=0,
iplotorbitals=0, 
ievaluate=0,
ekin0pp=0.0,
vxn0=0.0,
vyn0=0.0,
vzn0=-1.0,
eproj=0.0,
vpx=0.0,
vpy=0.0,
vpz=-1.0,
taccel=0.0,
trequest=100000,
timefrac=0.98,
ehom0=0.0,
ehomx=0.0,
ehomy=0.0,
ehomz=1.0,
ihome=0,
scatterelectronenergy=0.0,
scatterelectronw=1.0,
scatterelectronvxn=0.0,
scatterelectronvyn=0.0,
scatterelectronvzn=1.0,
scatterelectronx=0.0,
scatterelectrony=0.0,
scatterelectronz=0.0,
drcharges=5.0,
iflocaliz=0,                       
ismax=500,
itmax=10,
istinf=1,
ipasinf=1,
idyniter=0,        
iffastpropag=1,
ifexpevol=0,
irest=1,
istat=0, 
isave=1,
idenspl=0,
i3dz=0,
i3dx=0,
i3dstate=0,
modrho=999999,
jpos=0,
jvel=0,
jener=1,
jesc=0,
jforce=0,
jposcm=0,
jgeomion=0,
jinfo=1,
jdip=1,
jquad=0,
jang=0,
jspdp=0,
jenergy=1,
jgeomel=0,
jangabso=0,
jelf=0,
jstinf=1,
jstboostinv=0,
jstateoverlap=0,
nabsorb=0,
ifsicp=2,
ifredmas=0,
ionmdtyp=1,
icooltyp=0,
init_lcao=0,
ipsptyp=3,
ivdw=0,
izforcecorr=-1,
icooltimes=0,
ntref=0,
jheatmod=0,         
ifrhoint_time=0,
iangmo=0,
ifspemoms=0,
iftransme=0,
rheattemp=0.0,        
tempion=0.0,
dt1=0.001,
centfx=0.0,
centfy=0.0,
centfz=0.0,
shiftinix=0.0,
shiftiniy=0.0,
shiftiniz=0.0,
projvelx=0.0,
projvely=0.0,
projvelz=0.0,
projinix=0.0,
projiniy=0.0,
projiniz=0.0,
bcol1=0.0,
bcol23=0.0,
dbcol=0.1,
betacol=0.99,
chgcol=0.0,
ntheta=0,
nphi=0,
ispherabso=1,           
iangabso=0,
nangtheta=1,
nangphi=1,
ipes=0, 
angthetah=3.141592654,
angthetal=0.0,
angphil=0.0,
angphih=6.28,
mzforce=0,
myforce=0,
mxforce=0,       
nrare=0,
nfix=0,                  
idielec=0,                     
itft=3,
tnode=0.0,
deltat=0.0,
tpeak=0.0,
omega=0.0,
e0=0.0,
tfs=0.0,  
e1x=1.0,
e1y=0.0,
e1z=0.0,
e2x=0.0,
e2y=0.0,
e2z=0.0,
phi=0.0,
ilas=0,
projcharge=0.0,                  
iPotStatic=0,
write_for005=1,
modecalc=None,
calls=1,
nclust=0,nion=None
):
	self.modecalc = modecalc
        if self.modecalc is None:
	    self.modecalc = 'static'
        print("mode %s \n" % self.modecalc) 
        self.calls=calls
        #if (self.modecalc != 'static'):
        #      lines = open('pwtelemancall', 'r').readlines()
        #      self.calls = int(lines[0])
        #      print "Hello"
        #      print self.calls
        #self.name = label + '%04d' % self.calls
        #self.label = self.name
	self.label = label
        self.xc = xc
        self.iPotStatic=iPotStatic
        self.converged = False
        self.write_for005_file = write_for005
        self.for005g = {}
        self.for005 = {}
        self.for005d = {}
        self.for005s = {}
        self.e_fermi = None
	self.kxbox=kxbox
	self.kybox=kybox
	self.kzbox=kzbox
	self.kstate=kstate
	self.nspdw=nspdw
	self.dx=dx
	self.dy=dy
	self.dz=dz
	self.iforce=iforce
	self.ipseudo=ipseudo
	self.epswf=epswf
	self.e0dmp=e0dmp
	self.epsorc=epsorc
	self.b2occ=b2occ
	self.gamocc=gamocc
	self.deocc=deocc
	self.osfac=osfac
	self.temp=temp
	self.occmix=occmix
	self.epsoro=epsoro
	self.isurf=isurf
	self.nc=nc
	self.nk=nk
	self.endcon=endcon
	self.radjel=radjel
	self.surjel=surjel
	self.itback=itback
	self.bbeta=bbeta
	self.gamma=gamma
	self.beta4=beta4
	self.dpolx=dpolx
	self.dpoly=dpoly
	self.dpolz=dpolz
	self.iexcit=iexcit
	self.irotat=irotat
	self.ispidi=ispidi
	self.phirot=phirot
	self.scaleclust=scaleclust
	self.scaleclustx=scaleclustx
	self.scaleclusty=scaleclusty
	self.scaleclustz=scaleclustz
	self.iemomsrel=iemomsrel
	self.shiftclustx=shiftclustx
	self.shiftclusty=shiftclusty
	self.shiftclustz=shiftclustz
	self.rotclustx=rotclustx
	self.rotclusty=rotclusty
	self.rotclustz=rotclustz
	self.imob=imob
	self.shiftwfx=shiftwfx
	self.shiftwfy=shiftwfy
	self.shiftwfz=shiftwfz
	self.ispinsep=ispinsep
	self.phangle=phangle
	self.phphase=phphase
	self.idebug=idebug
	self.ifreezekspot=ifreezekspot
	self.itof=itof
	self.jescmask=jescmask
	self.ishiftcmtoorigin=ishiftcmtoorigin
	self.jescmaskorb=jescmaskorb
	self.ishutdown=ishutdown
	self.jplotdensitydiff=jplotdensitydiff
	self.jplotdensity2d=jplotdensity2d
	self.jplotdensitydiff2d=jplotdensitydiff2d
	self.nmptheta=nmptheta
	self.nmpphi=nmpphi
	self.jmp=jmp
	self.jovlp=jovlp
	self.jnorms=jnorms
	self.iscatterelectron=iscatterelectron
	self.jcharges=jcharges
	self.iaddcluster=iaddcluster
	self.iswforce=iswforce
	self.iplotorbitals=iplotorbitals
	self.ievaluate=ievaluate
	self.ekin0pp=ekin0pp
	self.vxn0=vxn0
	self.vyn0=vyn0
	self.vzn0=vzn0
	self.eproj=eproj
	self.vpx=vpx
	self.vpy=vpy
	self.vpz=vpz
	self.taccel=taccel
	self.trequest=trequest
	self.timefrac=timefrac
	self.ehom0=ehom0
	self.ehomx=ehomx
	self.ehomy=ehomy
	self.ehomz=ehomz
	self.ihome=ihome
	self.scatterelectronenergy=scatterelectronenergy
	self.scatterelectronw=scatterelectronw
	self.scatterelectronvxn=scatterelectronvxn
	self.scatterelectronvyn=scatterelectronvyn
	self.scatterelectronvzn=scatterelectronvzn
	self.scatterelectronx=scatterelectronx
	self.scatterelectrony=scatterelectrony
	self.scatterelectronz=scatterelectronz
	self.drcharges=drcharges
	self.iflocaliz=iflocaliz
	self.ismax=ismax
	self.itmax=itmax
	self.isave=isave
	self.irest=irest
	self.istat=istat
	self.istinf=istinf
	self.ipasinf=ipasinf
	self.idyniter=idyniter
	self.iffastpropag=iffastpropag
	self.ifexpevol=ifexpevol
	self.idenspl=idenspl
	self.i3dz=i3dz
	self.i3dx=i3dx
	self.i3dstate=i3dstate
	self.modrho=modrho
	self.jpos=jpos
	self.jvel=jvel
	self.jener=jener
	self.jesc=jesc
	self.jforce=jforce
	self.jposcm=jposcm
	self.jgeomion=jgeomion
	self.jinfo=jinfo
	self.jdip=jdip
	self.jquad=jquad
	self.jang=jang
	self.jspdp=jspdp
	self.jenergy=jenergy
	self.jgeomel=jgeomel
	self.jangabso=jangabso
	self.jelf=jelf
	self.jstinf=jstinf
	self.jstboostinv=jstboostinv
	self.jstateoverlap=jstateoverlap
	self.nabsorb=nabsorb
	self.ifsicp=ifsicp
	self.ifredmas=ifredmas
	self.ionmdtyp=ionmdtyp
	self.icooltyp=icooltyp
	self.init_lcao=init_lcao
	self.ipsptyp=ipsptyp
	self.ivdw=ivdw
	self.izforcecorr=izforcecorr
	self.icooltimes=icooltimes
	self.ntref=ntref
	self.jheatmod=jheatmod
	self.ifrhoint_time=ifrhoint_time
	self.iangmo=iangmo
	self.ifspemoms=ifspemoms
	self.iftransme=iftransme
	self.rheattemp=rheattemp
	self.tempion=tempion
	self.dt1=dt1
	self.centfx=centfx
	self.centfy=centfy
	self.centfz=centfz
	self.shiftinix=shiftinix
	self.shiftiniy=shiftiniy
	self.shiftiniz=shiftiniz
	self.bcol1=bcol1
	self.bcol23=bcol23
	self.dbcol=dbcol
	self.betacol=betacol
	self.chgcol=chgcol
	self.ntheta=ntheta
	self.nphi=nphi
	self.ispherabso=ispherabso
	self.iangabso=iangabso
	self.nangtheta=nangtheta
	self.nangphi=nangphi
	self.ipes=ipes
	self.angthetah=angthetah
	self.angthetal=angthetal
	self.angphil=angphil
	self.angphih=angphih
	self.mzforce=mzforce
	self.myforce=myforce
	self.mxforce=mxforce
	self.nrare=nrare
	self.nfix=nfix
	self.idielec=idielec
	self.itft=itft
	self.tnode=tnode
	self.deltat=deltat
	self.tpeak=tpeak
	self.omega=omega
	self.e0=e0
	self.tfs=tfs
	self.e1x=e1x
	self.e1y=e1y
	self.e1z=e1z
	self.e2x=e2x
	self.e2y=e2y
	self.e2z=e2z
	self.phi=phi
	self.ilas=ilas
	self.projcharge=projcharge
	self.projvelx=projvelx
	self.projvely=projvely
	self.projvelz=projvelz
	self.projinix=projinix
	self.projiniy=projiniy
	self.projiniz=projiniz

    def update(self, atoms):
	    if (self.modecalc is 'static'):
		self.irest=0
		self.isave=1
		self.itmax=10
		self.dt1=0.0
		print("irest %d \n" % self.irest)    
		print("itmax %d \n" % self.itmax)
		print("dt1 %d \n" % self.dt1)
		print("isave %d \n" % self.isave)
	    if (self.modecalc is 'dynamic'):
		self.irest=1
		self.isave=10
		self.itmax=10
		self.dt1=0.001
		if (self.calls  is 1):
			self.irest=0
		print("irest %d \n" % self.irest)    
		print("itmax %d  \n" % self.itmax)
		self.itmax=self.isave*self.calls
		print("isave %d  \n" % self.isave)
		print("calls %d  \n" % self.calls)
		print("itmax %d  \n" % self.itmax)
	    if (self.modecalc is 'plasmon'):
		self.irest=0
		self.isave=1000
		self.itmax=100000
		self.dt1=0.001
		if (self.calls  is 1):
			self.irest=0
		print("irest %d \n" % self.irest)    
		print("itmax %d  \n" % self.itmax)
		print("isave %d  \n" % self.isave)
		print("calls %d  \n" % self.calls)
		print("itmax %d  \n" % self.itmax)
	    if (not self.converged or
	    	len(self.numbers) != len(atoms) or
	    	(self.numbers != atoms.get_atomic_numbers()).any()):
		    	self.initialize(atoms)
		    	self.calculate(atoms)
	    elif ((self.positions != atoms.get_positions()).any()):
	    		self.calculate(atoms)
	    if (not self.converged or
	    	len(self.numbers) != len(atoms) or
	    	(self.numbers != atoms.get_atomic_numbers()).any()):
		    	self.initialize(atoms)
		    	self.calculate(atoms)
	    elif ((self.positions != atoms.get_positions()).any()):
	    		self.calculate(atoms)

    def initialize(self, atoms):
        self.numbers = atoms.get_atomic_numbers().copy()
        self.species = []
        for a, Z in enumerate(self.numbers):
            if Z not in self.species:
                self.species.append(Z)

        if 'pwteleman_PP_PATH' in os.environ:
            pppaths = os.environ['pwteleman_PP_PATH'].split(':')
        else:
            pppaths = []


        self.converged = False

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)

        if force_consistent:
            return self.etotal
        else:
            # Energy extrapolated to zero Kelvin:
            return  (self.etotal) 

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()


    def get_dipole_moment(self, atoms):
        return self.dipole

    def read_dipole(self):
        dipolemoment = np.zeros([1, 3])
        text = open('pdip.' + self.label, 'r').read().lower()
        for line in  iter(text.split('\n')):
            	if line.rfind('   0.00000 ') > -1:
                        dipolemoment = np.array([float(f) for f in line.split()[2:4]])
        return dipolemoment


    def calculate(self, atoms):
        self.positions = atoms.get_positions().copy()
        self.write_for005(atoms)

        pwteleman = os.environ['pwteleman_SCRIPT']
        locals = {'label': self.label}
        execfile(pwteleman, {}, locals)
        if (self.modecalc != 'static'):
              self.calls +=1
        fcalls=open('pwtelemancall','w')
        fcalls.write('%d \n' % self.calls)
        print('calls %d \n' % self.calls)
        fcalls.flush()
        fcalls.close()
        exitcode = locals['exitcode']
        if exitcode != 0:
        	execfile(pwteleman, {}, locals)
        	exitcode = locals['exitcode']
                if exitcode != 0:
       	       		 raise RuntimeError(('after retry  pwteleman exited with exit code: %d.  ' +
                                'Check %s.txt for more information.') %
                               (exitcode, self.label))

        self.dipole = self.read_dipole()
        self.read()
        self.converged = False

    def set_for005(self, key, value):
        """Set for005 parameter."""
        self.for005[key] = value

    def write_for005(self, atoms):
        """Write input parameters to for005-file."""
        #if(self.calls !=1):
        #      lines = open('pwtelemancall', 'r').readlines()
        #      self.calls = int(lines[0])
        #      print "Hello"
        #      print self.calls
        #      self.label = "pwteleman"
        #      self.name = self.label + '%04d' % self.calls
        #      self.label = self.name
        f2 = open('for005', 'w')
        #print self.label
        f2.write('%s \n' %self.label)
        f2.flush()
        f2.close()

        fh = open('for005.' +self.label, 'w')
        fion = open('for005ion.' +self.label, 'w')
        fdx = open('dx','w')
	fdx.write('%s \n' %self.dx)
	fdx.flush()
	fdx.close()
        fdx = open('nx','w')
	fdx.write('%s \n' %self.kxbox)
	fdx.flush()
	fdx.close()
        for005g = {
        'nclust':0,
	'nion':len(atoms),
'kxbox': self.kxbox,
'kybox': self.kybox,
'kzbox': self.kzbox,
'kstate': self.kstate,
'nspdw':self.nspdw,
'dx': self.dx,
'dy': self.dy,
'dz': self.dz,
'epswf': self.epswf,
'e0dmp': self.e0dmp,
'b2occ': self.b2occ,
'gamocc': self.gamocc,
'deocc': self.deocc,
'osfac': self.osfac,
'temp': self.temp,
'epsoro': self.epsoro,
'isurf': self.isurf,
'nc': self.nc,
'nk': self.nk,
'endcon': self.endcon,
'radjel': self.radjel,
'surjel': self.surjel,
'itback': self.itback,
'bbeta': self.bbeta,
'gamma': self.gamma,
'dpolx': self.dpolx,
'dpoly': self.dpoly,
'dpolz': self.dpolz,
'rotclustx': self.rotclustx,
'rotclusty': self.rotclusty,
'rotclustz': self.rotclustz,
'imob': self.imob,
'ispinsep': self.ispinsep,
'ishiftcmtoorigin': self.ishiftcmtoorigin,
'iaddcluster': self.iaddcluster,
'ievaluate': self.ievaluate,
'ehom0': self.ehom0,
'ehomx': self.ehomx,
'ehomy': self.ehomy,
'ehomz': self.ehomz,
'ihome': self.ihome,
'init_lcao': self.init_lcao,
'dbcol': self.dbcol,
'chgcol': self.chgcol,
'ntheta': self.ntheta,
'nphi': self.nphi,
	}
        for005d = {
'iforce': self.iforce,
'ipseudo': self.ipseudo,
'iexcit': self.iexcit,
'irotat': self.irotat,
'ispidi': self.ispidi,
'phirot': self.phirot,
'iemomsrel': self.iemomsrel,
'phangle': self.phangle,
'phphase': self.phphase,
'ifreezekspot': self.ifreezekspot,
'itof': self.itof,
'jescmask': self.jescmask,
'jescmaskorb': self.jescmaskorb,
'jplotdensitydiff': self.jplotdensitydiff,
'jplotdensity2d': self.jplotdensity2d,
'jplotdensitydiff2d': self.jplotdensitydiff2d,
'nmpphi': self.nmpphi,
'jmp': self.jmp,
'jovlp': self.jovlp,
'jnorms': self.jnorms,
'iscatterelectron': self.iscatterelectron,
'jcharges': self.jcharges,
'ekin0pp': self.ekin0pp,
'vxn0': self.vxn0,
'vyn0': self.vyn0,
'vzn0': self.vzn0,
'eproj': self.eproj,
'vpx': self.vpx,
'vpy': self.vpy,
'vpz': self.vpz,
'taccel': self.taccel,
'trequest': self.trequest,
'timefrac': self.timefrac,
'scatterelectronenergy': self.scatterelectronenergy,
'scatterelectronw': self.scatterelectronw,
'scatterelectronvxn': self.scatterelectronvxn,
'scatterelectronvyn': self.scatterelectronvyn,
'scatterelectronvzn': self.scatterelectronvzn,
'scatterelectronx': self.scatterelectronx,
'scatterelectrony': self.scatterelectrony,
'scatterelectronz': self.scatterelectronz,
'drcharges': self.drcharges,
'iflocaliz': self.iflocaliz,
'ismax': self.ismax,
'itmax': self.itmax,
'istinf': self.istinf,
'ipasinf': self.ipasinf,
'idyniter': self.idyniter,
'iffastpropag': self.iffastpropag,
'ifexpevol': self.ifexpevol,
'irest': self.irest,
'istat': self.istat,
'isave': self.isave,
'idenspl': self.idenspl,
'i3dz': self.i3dz,
'i3dx': self.i3dx,
'i3dstate': self.i3dstate,
'modrho': self.modrho,
'jpos': self.jpos,
'jvel': self.jvel,
'jener': self.jener,
'jesc': self.jesc,
'jforce': self.jforce,
'jposcm': self.jposcm,
'jgeomion': self.jgeomion,
'jinfo': self.jinfo,
'jdip': self.jdip,
'jquad': self.jquad,
'jang': self.jang,
'jspdp': self.jspdp,
'jenergy': self.jenergy,
'jgeomel': self.jgeomel,
'jelf': self.jelf,
'jstinf': self.jstinf,
'jstboostinv': self.jstboostinv,
'nabsorb': self.nabsorb,
'ifsicp': self.ifsicp,
'ifredmas': self.ifredmas,
'ionmdtyp': self.ionmdtyp,
'icooltyp': self.icooltyp,
'ipsptyp': self.ipsptyp,
'ivdw': self.ivdw,
'izforcecorr': self.izforcecorr,
'ntref': self.ntref,
'ifrhoint_time': self.ifrhoint_time,
'iangmo': self.iangmo,
'ifspemoms': self.ifspemoms,
'iftransme': self.iftransme,
'tempion': self.tempion,
'dt1': self.dt1,
'centfx': self.centfx,
'centfy': self.centfy,
'centfz': self.centfz,
'shiftinix': self.shiftinix,
'shiftiniy': self.shiftiniy,
'shiftiniz': self.shiftiniz,
'ispherabso': self.ispherabso,
'iangabso': self.iangabso,
'nangtheta': self.nangtheta,
'nangphi': self.nangphi,
'ipes': self.ipes,
'angthetah': self.angthetah,
'angthetal': self.angthetal,
'angphil': self.angphil,
'angphih': self.angphih,
'mzforce': self.mzforce,
'myforce': self.myforce,
'mxforce': self.mxforce,
'nfix': self.nfix,
'itft': self.itft,
'tnode': self.tnode,
'deltat': self.deltat,
'tpeak': self.tpeak,
'omega': self.omega,
'e0': self.e0,
'e1x': self.e1x,
'e1y': self.e1y,
'e1z': self.e1z,
'e2x': self.e2x,
'e2y': self.e2y,
'e2z': self.e2z,
'phi': self.phi,
'projcharge': self.projcharge,
'projvelx': self.projvelx,
'projvely': self.projvely,
'projvelz': self.projvelz,
'projinix': self.projinix,
'projiniy': self.projiniy,
'projiniz': self.projiniz,
        'ifhamdiag':0
            }
        for005s = {
}

        for005g.update(self.for005g)
        for005d.update(self.for005d)
        for005s.update(self.for005s)

        fh.write('&GLOBAL\n')
        for key, value in for005g.items():
            if value is None:
                continue

            if isinstance(value, list):
                for line in value:
                    fh.write(line + ',\n')
            else:
                unit = keys_with_units.get(for005ify(key))
                if unit is None:
                    fh.write('%s=%s,\n' % (key, value))
                else:
                    if 'fs**2' in unit:
                        value /= fs**2
                    elif 'fs' in unit:
                        value /= fs
                    fh.write('%s=%f,\n' % (key, value))
        fh.write('&END\n')
        fh.write('&DYNAMIC\n')
        for key, value in for005d.items():
            if value is None:
                continue

            if isinstance(value, list):
                for line in value:
                    fh.write(line + ',\n')
            else:
                unit = keys_with_units.get(for005ify(key))
                if unit is None:
                    fh.write('%s=%s,\n' % (key, value))
                else:
                    if 'fs**2' in unit:
                        value /= fs**2
                    elif 'fs' in unit:
                        value /= fs
                    fh.write('%s=%f,\n' % (key, value))
        fh.write('&END\n')
        fh.write('&SURFACE\n')
        for key, value in for005s.items():
            if value is None:
                continue

            if isinstance(value, list):
                for line in value:
                    fh.write(line + ',\n')
            else:
                unit = keys_with_units.get(for005ify(key))
                if unit is None:
                    fh.write('%s=%s,\n' % (key, value))
                else:
                    if 'fs**2' in unit:
                        value /= fs**2
                    fh.write('%s=%f,\n' % (key, value))
        fh.write('&END\n')

        a = 0
        for pos, Z in zip(self.positions, self.numbers):
            pos = pos / 0.52
            spin = a % 2 *2-1
            a += 1
#            fion.write('%.14f %.14f %.14f 1 xyz 1,., ' %  tuple(pos))
#            fion.write('%i\n' %  spin)
            fion.write('%.14f %.14f %.14f ' %  tuple(pos))
            fion.write('%i xyz  1.0 ' %  Z)
            fion.write('%i\n' %  spin)

        fh.flush()
        fion.flush()
        fh.close()
        fion.close()



    def read(self):
        """Read results from pwteleman's text-output file."""
        text = open('energies.' +self.label, 'r').readlines()
        while (len(text)==0):
          text = open('energies.' +self.label, 'r').readlines()

        for line in  text:
                        self.etotal = float(line.split('=')[-1])*Hartree
                        break
        # Forces (changed so forces smaller than -999eV/A can be fetched):
        lines= open('forces.' +self.label, 'r').readlines()
        while (len(lines)!=3*len(self.numbers)):
          lines= open('forces.' +self.label, 'r').readlines()

        self.forces = np.zeros((len(self.numbers), 3))
        for i in range(len(self.numbers)):
            self.forces[i, 0] = float(lines[3*i])
            self.forces[i, 1] = float(lines[3*i+1])
            self.forces[i, 2] = float(lines[3*i+2])



def getrecord(fileobj, dtype):
    """Used to read in binary files.
    """
    typetosize = {'l':4, 'f':4, 'd':8}# XXX np.int, np.float32, np.float64
    assert dtype in typetosize # XXX
    size = typetosize[dtype]
    record = array.array('l')
    trunk = array.array(dtype)
    record.fromfile(fileobj, 1)
    nofelements = int(record[-1]) / size
    trunk.fromfile(fileobj, nofelements)
    record.fromfile(fileobj, 1)
    data = np.array(trunk, dtype=dtype)
    if len(data)==1:
        data = data[0]
    return data



def for005ify(key):
    return key.lower().replace('_', '').replace('.', '').replace('-', '')


keys_with_units = {
	'temp':'Hartree',
	'dx':'Bohr',
	'dy':'Bohr',
	'dz':'Bohr',
	'b2occ':'Hartree',
	'deocc':'Hartree',
	'epsoro':'Hartree',
	'epswf':'Hartree',
	'e0dmp':'Hartree',
        'dt1': 'fs',
        'e0':'Hartree'}
    
