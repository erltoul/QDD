'''
This file is a part of PW-TELEMAN project.
PW-TELEMAN is a Time-Dependent Electronic Dynamics in Molecules And Nanosystems library.
Copyright (C) 2011-2015  Paul-Gerhard Reinhard, Eric Suraud, Florent Calvayrac,
Phuong Mai Dinh, David Brusson, Philipp Wopperer, José María Escartín Esteban.

PW-Teleman is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PW-Teleman is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PW-Teleman.  If not, see <http://www.gnu.org/licenses/>.
'''

import os
import time
pwteleman = '$TELEMANROOT/code/telemanexec'
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
