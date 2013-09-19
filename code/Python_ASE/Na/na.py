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

na13=BodyCenteredCubic('Na', [(1, 0, 0)], [1], latticeconstant=3.8 * sqrt(2))
write('na13.xyz',na13)
#from ase.calculators.pwteleman import pwteleman
from ase.calculators.pwtelemandynr import pwtelemandynr

na13.set_calculator(pwtelemandynr(isurf=0,modecalc='static'))
#e = na13.get_potential_energy()
#print e
traj=PickleTrajectory('na13.traj','w')
#dyn=QuasiNewton(na13).run(fmax=0.01)
e = na13.get_potential_energy()
print e
na13.set_calculator(pwtelemandynr(isurf=0,modecalc='dynamic'))
dyn = FIRE(na13, dt=0.1)
dyn.run(fmax=0.01,steps=3)


