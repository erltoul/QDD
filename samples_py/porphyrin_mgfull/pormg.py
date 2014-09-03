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

pormg=read('pormg.xyz')

pormg.set_calculator(pwtelemandynr(kxbox=96,kybox=96,kzbox=96,init_lcao=1,ismax=100,modecalc='static'))
e = pormg.get_potential_energy()
traj=PickleTrajectory('pormg.traj','w')
dyn=QuasiNewton(pormg).run(fmax=0.5)
write('newpormg.xyz',pormg)
pormg.set_calculator(pwtelemandynr(kxbox=96,kybox=96,kzbox=96,init_lcao=1,ismax=100,itmax=148809,dt1=0.001,centfx=0.1,centfy=0.1,centfz=0.1,modecalc='plasmon'))
e = pormg.get_potential_energy()
