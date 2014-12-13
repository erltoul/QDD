import subprocess
import ase.utils.geometry as geometry
import ase.io as io
from ase.lattice.spacegroup import crystal
from ase.visualize import view
from ase.io import write
from ase.lattice.surface import surface
from ase.data.molecules import molecule
from ase.lattice.surface import add_adsorbate
from ase.io import read


from math import sqrt
from ase.cluster.cubic import BodyCenteredCubic
from ase.tasks.molecule import MoleculeTask
from ase.data import covalent_radii, atomic_numbers

from ase import Atoms, Atom
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton,FIRE
from ase.io import PickleTrajectory
from ase.neb import NEB
from ase.calculators.emt import EMT

from pwtelemandynr import pwtelemandynr


na13=BodyCenteredCubic('Na', [(1, 0, 0)], [1], latticeconstant=3.8 * sqrt(2))
write('na13.xyz',na13)
box=128
dx0=0.8
l=3
for m in range(0,l):
		print m
		print "setup"
		rlabel='%04d' % m
		print rlabel
		bbox=box/pow(2,m)
		dx1=dx0*pow(2,m) 
		na13.set_calculator(pwtelemandynr(init_lcao=0,label=rlabel,ismax=15,kxbox=bbox,kybox=bbox,kzbox=bbox,dx=dx1,dy=dx1,dz=dx1,modecalc='static',ipsptyp=0))
		e1hp = na13.get_potential_energy()
		print e1hp


print "end of setup"

for ev in range(1,4):
	for m in range(1,l):
		print m
		print "coarser"
		rlabel='%04d' % m
		rlabelm='%04d' % (m-1)
		bbox=box/pow(2,m)
		dx1=dx0*pow(2,m) 
	        subprocess.call(['cp',rlabelm+'/rsave.'+rlabelm,rlabel+'/rsave.'+rlabel],shell=False)	
	        print rlabelm+'/rsave.'+rlabelm
	        print rlabel+'/rsave.'+rlabel
		na13.set_calculator(pwtelemandynr(istat=1,init_lcao=0,label=rlabel,ismax=40*ev,kxbox=bbox,kybox=bbox,kzbox=bbox,dx=dx1,dy=dx1,dz=dx1,modecalc='static_coarse',ipsptyp=0))
		e1hp = na13.get_potential_energy()
		print e1hp
	for m in range(l-1,-1,-1):
		print m
		print "finer"
		rlabel='%04d' % m
		rlabelm='%04d' % (m+1)
		bbox=box/pow(2,m)
		dx1=dx0*pow(2,m) 
	        subprocess.call(['cp',rlabelm+'/rsave.'+rlabelm,rlabel+'/rsave.'+rlabel],shell=False)	
	        print rlabelm+'/rsave.'+rlabelm
	        print rlabel+'/rsave.'+rlabel
		na13.set_calculator(pwtelemandynr(istat=1,init_lcao=0,label=rlabel,ismax=40*ev,kxbox=bbox,kybox=bbox,kzbox=bbox,dx=dx1,dy=dx1,dz=dx1,modecalc='static_fine',ipsptyp=0))
		e1hp = na13.get_potential_energy()
		print m
		print e1hp


