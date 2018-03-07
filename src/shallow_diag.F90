!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Copyright 2007.  Los Alamos National Security, LLC. This material was
!produced under U.S. Government contract DE-AC52-06NA25396 for Los
!Alamos National Laboratory (LANL), which is operated by Los Alamos
!National Security, LLC for the U.S. Department of Energy. The
!U.S. Government has rights to use, reproduce, and distribute this
!software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
!LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
!FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
!derivative works, such modified software should be clearly marked, so
!as not to confuse it with the version available from LANL.
!
!Additionally, this program is free software; you can redistribute it
!and/or modify it under the terms of the GNU General Public License as
!published by the Free Software Foundation; either version 2 of the
!License, or (at your option) any later version. Accordingly, this
!program is distributed in the hope that it will be useful, but WITHOUT
!ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
!for more details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "macros.h"

!
!  diagnostics for shallow water model
!
subroutine output_model(doit_model,doit_diag,time,Q,Qhat,q1,q2,q3,work1,work2)
use params
use spectrum
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(*)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)
real*8 :: q3(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time
logical :: doit_model,doit_diag
integer i,j,k,n,ierr,csig
integer :: n1,n1d,n2,n2d,n3,n3d,pv_type,stype

real*8 :: tsave

character(len=80) :: message
character,save :: pow_access="0"
CPOINTER fid,fidj,fidS,fidcore

! append to output files, unless this is first call.                                                                                       
if (pow_access=="0") then
   pow_access="w"
else
   pow_access="a"
endif

if (compute_transfer) then
   ! spec_r computed last time step
   ! spec_diff,spec_model were computed in RHS computation at the
   ! beginning of this flag (becuase compute_transfer flag was set)
   ! So they are all known at time_old. Now compute spec_r_new 
   ! (used to compute edot_r at mid-time level)
   call compute_Edotspec_shallow(time,Q,q1,work1,work2)

   ! e_dot_r    known at (time + time_old)/2
   ! spec_diff,spec_model  known at time_old
   ! compute spec_diff at 'time' be re-calling getrhs with compute_transfer
   ! flag still set = .true. and compute_ints==1

   spec_diff_new=spec_diff  
   spec_model_new=spec_model
   tsave=transfer_comp_time
   call getrhs(4,q1,Qhat,Q,time,1,work1,work2)
   transfer_comp_time=tsave

   spec_diff = (spec_diff_new+spec_diff)/2
   spec_model = (spec_model_new+spec_model)/2
   
   ! output all the spectrum:
   call output_tran(time,Q,q1,q2,q3,work1,work2)


   compute_transfer=.false.
endif


if (doit_model .or. time==time_initial) then
if ( g_bdy_x1==PERIODIC .and. &
     g_bdy_y1==PERIODIC .and. &
     g_bdy_z1==PERIODIC) then
   call compute_spec_shallow(time,Q,q1,work1,work2)
   call output_spec(time,time_initial)
   
   !set this flag so that for next timestep, we will compute and save
   !spectral transfer functions:
   compute_transfer=.true.
endif
endif


!if (doit_diag == .true.) then ! If you want it every time step comment this out.
! output data for time series
if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time_initial
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".kepower"
   call copen(message,pow_access,fid,ierr)
!   write(6,*) "Opening power spectrum file"
   if (ierr/=0) then
      write(message,'(a,i5)') "power_output(): Error opening .kepower file errno=",ierr
      call abortdns(message)
   endif
   call cwrite8(fid,time, 1)
   call cwrite8(fid,Q(5,5,1,1), 1)
   call cwrite8(fid,Q(5,5,1,2), 1)
   call cwrite8(fid,Q(5,5,1,3), 1)
   call cclose(fid,ierr)
endif
!endif

end subroutine
