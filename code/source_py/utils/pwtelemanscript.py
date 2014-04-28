import os
import time
pwteleman = 'mpirun -np 2 ~/gitcalv/code/essai.par'
#pwteleman = '~/gitcalv/code/essai.seq'
line = open('for005', 'r').readlines()
label = line[0].split(" \n")

#exitcode=os.system(' touch lock.'+label[0])
exitcode=os.system(' touch pwtelemancall')
exitcode=os.system(' touch pdip.'+label[0])
exitcode=os.system(' touch energies.'+label[0])
exitcode=os.system(' touch forces.'+label[0])
try:
  line = open('lock.' +label[0], 'r').readlines()
  test=0
except IOError:
  test=1

#print test
time.sleep(0)
if(test==0):
  while (line[0]=="lock"):
    line = open('lock.' +label[0], 'r').readlines()
exitcode=os.system(' %s  > outputpwteleman ' % (pwteleman))
