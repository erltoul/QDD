!This file is a part of PW-TELEMAN project.
!PW-TELEMAN is a Time-Dependent Electronic Dynamics in Molecules And Nanosystems library.
!Copyright (C) 2011-2015  Paul-Gerhard Reinhard, Eric Suraud, Florent Calvayrac,
!Phuong Mai Dinh, David Brusson, Philipp Wopperer, José María Escartín Esteban.
!
!PW-Teleman is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!PW-Teleman is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with PW-Teleman.  If not, see <http://www.gnu.org/licenses/>.

PROGRAM Create_wisdom

!This program compute optimal FFTW plans for some usual grid sizes
!For 1d plans, set def_1d to 1 in define.h
!For 3d plans, set def_3d to 1 in define.h
!Both oh them can be set to 1d
!For a specific plan, set spec_plan to 1 in define.h, then define if it's a 3d or 1d plan (or both) 
!with def_1d or 3d flags. Finally, set the size with nxspec, nyspec and nzspec in this program (only nxspec for 1d)

#include"define.h"
USE, intrinsic :: iso_c_binding
USE FFTW

type(C_PTR) :: wisdom_plan1df,wisdom_plan3df,wisdom_plan1db,wisdom_plan3db
INTEGER(C_INT) :: wisdomtest
COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: mat1d(:),mat3d(:,:,:) ! 2d is not used in the tddft code
INTEGER :: grid_size1d(7)=(/24,48,64,96,128,192,256/),grid_size3d(5)=(/24,48,64,96,128/)
NAMELIST /gridsize/ nxspec,nyspec,nzspec

wisdomtest=fftw_import_wisdom_from_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
IF (wisdomtest == 0) WRITE(6,*) 'wisdom_fftw.dat not found, creating it'

#if(spec_plan)
OPEN(UNIT=16,STATUS='old',FORM='formatted',FILE='spec_plan.dat')
READ(16,gridsize)
CLOSE(16)
WRITE(6,*)'GRIDSIZE READ'
WRITE(6,*)'nxspec,nyspec,nzspec =',nxspec,nyspec,nzspec
#if(def_1d)
ALLOCATE(mat1d(nxspec))
wisdom_plan1df=fftw_plan_dft_1d(nxspec,mat1d,mat1d,FFTW_FORWARD,FFTW_EXHAUSTIVE)
wisdom_plan1db=fftw_plan_dft_1d(nxspec,mat1d,mat1d,FFTW_BACKWARD,FFTW_EXHAUSTIVE)
wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
IF (wisdomtest == 0)   WRITE(6,*) 'Error exporting wisdom to file wisdom_fftw.dat'
DEALLOCATE(mat1d)
CALL fftw_destroy_plan(wisdom_plan1df)
CALL fftw_destroy_plan(wisdom_plan3df)
WRITE(6,*)'Plan 1d',nxspec,'done'
#endif !def_1d

#if(def_3d)
ALLOCATE(mat3d(nxspec,nyspec,nzspec))
wisdom_plan3df=fftw_plan_dft_3d(nxspec,nyspec,nzspec,mat3d,mat3d,FFTW_FORWARD,FFTW_EXHAUSTIVE)
wisdom_plan3db=fftw_plan_dft_3d(nxspec,nyspec,nzspec,mat3d,mat3d,FFTW_BACKWARD,FFTW_EXHAUSTIVE)
wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
IF (wisdomtest == 0)   WRITE(6,*) 'Error exporting wisdom to file wisdom_fftw.dat'
DEALLOCATE(mat3d)
CALL fftw_destroy_plan(wisdom_plan3df)
CALL fftw_destroy_plan(wisdom_plan3db)
WRITE(6,*)'Plan 3d',nxspec,nyspec,nzspec,'done'
#endif !def_3d

#else

#if(def_1d)
DO ii=1,size(grid_size1d)
  nb=grid_size1d(ii)
  ALLOCATE(mat1d(nb))
  wisdom_plan1df=fftw_plan_dft_1d(nb,mat1d,mat1d,FFTW_FORWARD,FFTW_EXHAUSTIVE)
  wisdom_plan1db=fftw_plan_dft_1d(nb,mat1d,mat1d,FFTW_BACKWARD,FFTW_EXHAUSTIVE)
  wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
  IF (wisdomtest == 0)   WRITE(6,*) 'Error exporting wisdom to file wisdom_fftw.dat'
  DEALLOCATE(mat1d)
  CALL fftw_destroy_plan(wisdom_plan1df)
  CALL fftw_destroy_plan(wisdom_plan1db)
  WRITE(6,*)'Plan 1d',nb,'done'
END DO
#endif !def_1d

#if(def_3d)
DO ii=1,size(grid_size3d)
  nb=grid_size3d(ii)
  ALLOCATE(mat3d(nb,nb,nb))
  wisdom_plan3df=fftw_plan_dft_3d(nb,nb,nb,mat3d,mat3d,FFTW_FORWARD,FFTW_EXHAUSTIVE)
  wisdom_plan3db=fftw_plan_dft_3d(nb,nb,nb,mat3d,mat3d,FFTW_BACKWARD,FFTW_EXHAUSTIVE)
  wisdomtest=fftw_export_wisdom_to_filename(C_CHAR_'wisdom_fftw.dat'//C_NULL_CHAR)
  IF (wisdomtest == 0)   WRITE(6,*) 'Error exporting wisdom to file wisdom_fftw.dat'
  DEALLOCATE(mat3d)
  CALL fftw_destroy_plan(wisdom_plan3df)
  CALL fftw_destroy_plan(wisdom_plan3db)
  WRITE(6,*)'Plan 3d',nb,nb,nb,'done'
END DO
#endif !def_3d
#endif !spec_plan

END PROGRAM Create_wisdom
