!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.

!*******************************************************************************
module io
!*******************************************************************************
! Navigation map:
!   - openfiles: output/restart file setup
!   - energy and wall-stress writers: lightweight diagnostics
!   - CGNS writers: optional PPCGNS output paths
!   - output_loop and inst_write: periodic field/statistics output
!   - force_tot: integrated force diagnostics
!   - checkpoint/output_final/output_init: restart and final output
!
! Keep binary format changes documented in data-ingestion or validation notes;
! downstream surrogate/export tooling depends on these output conventions.
use types, only : rprec
use param, only : ld, nx, ny, nz, nz_tot, path, coord, rank, nproc, jt_total
use param, only : total_time, total_time_dim, lbz, jzmin, jzmax
use param, only : cumulative_time
use sim_param , only : w, dudz, dvdz
use sgs_param , only : Cs_opt2
use string_util, only : string_concat, string_splice
use messages, only : error, warn
use time_average, only : tavg_t
#ifdef PPMPI
use mpi
#endif

#ifdef PPCGNS
use cgns
#ifdef PPMPI
use param, only: ierr
#endif
#endif

implicit none
save
private

public jt_total, openfiles, energy, output_loop, output_final, output_init,    &
    write_tau_wall_bot, write_tau_wall_top

! Where to end with nz index.
integer :: nz_end

! time averaging
type(tavg_t) :: tavg

character(:), allocatable :: fcumulative_time

contains

!*******************************************************************************
subroutine openfiles()
!*******************************************************************************
use param, only : use_cfl_dt, dt, cfl_f, checkpoint_file
use sgs_param, only : lagran_dt
implicit none
logical :: exst
integer :: ios
character(len=512) :: restart_line

! Time-step and CFL values read from cumulative-time restart metadata.
real(rprec) :: dt_r, cfl_r, lagran_dt_r

! Create file names
allocate(fcumulative_time, source = path // 'total_time.dat')
allocate(checkpoint_file , source = path // 'vel.out')
dt_r = dt
cfl_r = cfl_f
lagran_dt_r = 0._rprec

if (cumulative_time) then
    inquire (file=fcumulative_time, exist=exst)
    if (exst) then
        open (1, file=fcumulative_time)
        read(1, '(a)', iostat=ios) restart_line
        close (1)
        if (ios /= 0) then
            call error('openfiles', 'Could not read cumulative restart metadata')
        else
            read(restart_line, *, iostat=ios) jt_total, total_time,             &
                total_time_dim, dt_r, cfl_r, lagran_dt_r
            if (ios /= 0) then
                ! Backward compatibility with the historical five-value file.
                read(restart_line, *) jt_total, total_time, total_time_dim,     &
                    dt_r, cfl_r
                lagran_dt_r = 0._rprec
                if (coord == 0) call warn('openfiles',                         &
                    'Legacy total_time.dat has no lagran_dt SGS history; ' //  &
                    'using zero until the next dynamic update')
            end if
            lagran_dt = lagran_dt_r
        end if
    else
        ! assume this is the first run on cumulative time
        if ( coord == 0 ) then
            write (*, *) '--> Assuming jt_total = 0, total_time = 0.0'
        end if
        jt_total = 0
        total_time = 0._rprec
        total_time_dim = 0._rprec
    end if
end if

! Update dynamic time stepping info if required; otherwise discard.
if ( use_cfl_dt ) then
    dt = dt_r
    cfl_f = cfl_r
end if

end subroutine openfiles

!*******************************************************************************
subroutine energy (ke)
!*******************************************************************************
use types, only : rprec
#ifdef PPMPI
use param, only : comm, ierr, MPI_RPREC
#endif
use sim_param, only : u, v, w
implicit none
integer :: jx, jy, jz, nan_count
real(rprec)::KE
#ifdef PPMPI
real(rprec) :: ke_global
#endif

! Initialize variables
nan_count = 0
ke = 0._rprec

#if defined(PPLES_GPU)
! Device-side reduction over the resident u,v,w - replaces the full-field
! D2H + host loop. Same summand; only the (diagnostic-only) parallel
! reduction order differs from the serial host sum.
!$acc parallel loop collapse(3) default(present) reduction(+:ke)
do jz = 1, nz-1
do jy = 1, ny
do jx = 1, nx
    ke = ke + u(jx,jy,jz)*u(jx,jy,jz)                                         &
        + v(jx,jy,jz)*v(jx,jy,jz)                                             &
        + 0.25_rprec*(w(jx,jy,jz)+w(jx,jy,jz+1))                              &
        * (w(jx,jy,jz)+w(jx,jy,jz+1))
end do
end do
end do
#else
do jz = 1, nz-1
do jy = 1, ny
do jx = 1, nx
    ke = ke + u(jx,jy,jz)*u(jx,jy,jz)                                         &
        + v(jx,jy,jz)*v(jx,jy,jz)                                             &
        + 0.25_rprec*(w(jx,jy,jz)+w(jx,jy,jz+1))                              &
        * (w(jx,jy,jz)+w(jx,jy,jz+1))
end do
end do
end do
#endif

! Perform spatial averaging
ke = ke*0.5_rprec/(nx*ny*(nz-1))

#ifdef PPMPI
call mpi_reduce (ke, ke_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
if (rank == 0) then  ! note that it's rank here, not coord
    ke = ke_global/nproc
#endif
    open(2,file=path // 'output/check_ke.dat', status='unknown',               &
        form='formatted', position='append')
    write(2,*) total_time,ke
    close(2)
#ifdef PPMPI
end if
#endif

end subroutine energy

!*******************************************************************************
subroutine write_tau_wall_bot()
!*******************************************************************************
use types ,only: rprec
use param ,only: jt_total, total_time, total_time_dim, dt, dt_dim, wbase
use param ,only: L_x, z_i, u_star
use functions ,only: get_tau_wall_bot
implicit none

real(rprec) :: turnovers

turnovers = total_time_dim / (L_x * z_i / u_star)

open(2,file=path // 'output/tau_wall_bot.dat', status='unknown',               &
    form='formatted', position='append')

!! one time header output
if (jt_total==wbase) write(2,*)                                                &
    'jt_total, total_time, total_time_dim, turnovers, dt, dt_dim, 1.0, tau_wall'

!! continual time-related output
write(2,*) jt_total, total_time, total_time_dim, turnovers, dt, dt_dim,        &
    1.0, get_tau_wall_bot()
close(2)

end subroutine write_tau_wall_bot

!*******************************************************************************
subroutine write_tau_wall_top()
!*******************************************************************************
use types, only : rprec
use param, only : jt_total, total_time, total_time_dim, dt, dt_dim, wbase
use param, only : L_x, z_i, u_star
use functions, only : get_tau_wall_top
implicit none

real(rprec) :: turnovers

turnovers = total_time_dim / (L_x * z_i / u_star)

open(2,file=path // 'output/tau_wall_top.dat', status='unknown',               &
    form='formatted', position='append')

! one time header output
if (jt_total==wbase) write(2,*)                                                &
    'jt_total, total_time, total_time_dim, turnovers, dt, dt_dim, 1.0, tau_wall'

! continual time-related output
write(2,*) jt_total, total_time, total_time_dim, turnovers, dt, dt_dim,        &
    1.0, get_tau_wall_top()
close(2)

end subroutine write_tau_wall_top

#ifdef PPCGNS
#ifdef PPMPI
!*******************************************************************************
subroutine write_parallel_cgns (file_name, nx, ny, nz, nz_tot, start_n_in,     &
    end_n_in, xin, yin, zin, num_fields, fieldNames, input )
!*******************************************************************************
implicit none

integer, intent(in) :: nx, ny, nz, nz_tot, num_fields
! Name of file to be written
character(*), intent(in) :: file_name
! Name of fields we are writing
character(*), intent(in), dimension(:) :: fieldNames
! Data to be written
real(rprec), intent(in), dimension(:) :: input
! Coordinates to write
real(rprec), intent(in), dimension(:) :: xin, yin, zin
! Where the total node counter starts nodes
integer, intent(in) :: start_n_in(3)
! Where the total node counter ends nodes
integer, intent(in) :: end_n_in(3)

integer :: fn=1        ! CGNS file index number
integer :: ier         ! CGNS error status
integer :: base=1      ! base number
integer :: zone=1      ! zone number
integer :: nnodes      ! Number of nodes in this processor
integer :: sol =1      ! solution number
integer :: field       ! section number
integer(cgsize_t) :: sizes(3,3)  ! Sizes

! Convert input to right data type
integer(cgsize_t) :: start_n(3)  ! Where the total node counter starts nodes
integer(cgsize_t) :: end_n(3)  ! Where the total node counter ends nodes

! Building the lcoal mesh
integer :: i,j,k
real(rprec), dimension(nx,ny,nz) :: xyz

!  ! Set the parallel communicator
!  call cgp_mpi_comm_f(cgnsParallelComm, ierr)

! Convert types such that CGNS libraries can handle the input
start_n(1) = int(start_n_in(1), cgsize_t)
start_n(2) = int(start_n_in(2), cgsize_t)
start_n(3) = int(start_n_in(3), cgsize_t)
end_n(1) = int(end_n_in(1), cgsize_t)
end_n(2) = int(end_n_in(2), cgsize_t)
end_n(3) = int(end_n_in(3), cgsize_t)

! The total number of nodes in this processor
nnodes = nx*ny*nz

! Sizes, used to create zone
sizes(:,1) = (/int(nx, cgsize_t),int(ny, cgsize_t),int(nz_tot, cgsize_t)/)
sizes(:,2) = (/int(nx-1, cgsize_t),int(ny-1, cgsize_t),int(nz_tot-1, cgsize_t)/)
sizes(:,3) = (/int(0, cgsize_t) , int(0, cgsize_t), int(0, cgsize_t)/)

! Open CGNS file
call cgp_open_f(file_name, CG_MODE_WRITE, fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write base
call cg_base_write_f(fn, 'Base', 3, 3, base, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write zone
call cg_zone_write_f(fn, base, 'Zone', sizes, Structured, zone, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write print info to screen
if (coord .eq. 0) then
    write(*,*) 'Writing, ', file_name
end if

! Create data nodes for coordinates
call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateX', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateY', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateZ', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

! This is done for the 3 dimensions x,y and z
! It writes the coordinates
! Create grid points
do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = xin(i)
end do
end do
end do

call cgp_coord_write_data_f(fn, base, zone, 1,                                 &
    start_n, end_n, xyz(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = yin(j)
end do
end do
end do
call cgp_coord_write_data_f(fn, base, zone, 2,   &
    start_n, end_n, xyz(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = zin(k)
end do
end do
end do
call cgp_coord_write_data_f(fn, base, zone, 3,   &
                            start_n, end_n, xyz(1:nx,1:ny,1:nz), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Create a centered solution
call cg_sol_write_f(fn, base, zone, 'Solution', Vertex, sol, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the solution
do i=1,num_fields
    call cgp_field_write_f(fn, base, zone, sol, RealDouble, fieldNames(i),     &
        field, ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

    call cgp_field_write_data_f(fn, base, zone, sol, field, start_n, end_n,    &
        input((i-1)*nnodes+1:(i)*nnodes), ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

end do

! Close the file
call cgp_close_f(fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

end subroutine write_parallel_cgns

!*******************************************************************************
subroutine write_null_cgns (file_name, nx, ny, nz, nz_tot, start_n_in,         &
    end_n_in, xin, yin, zin, num_fields, fieldNames )
!*******************************************************************************
implicit none

integer, intent(in) :: nx, ny, nz, nz_tot, num_fields
! Name of file to be written
character(*), intent(in) :: file_name
! Name of fields we are writing
character(*), intent(in), dimension(:) :: fieldNames
! Coordinates to write
real(rprec), intent(in), dimension(:) :: xin, yin, zin
! Where the total node counter starts nodes
integer, intent(in) :: start_n_in(3)
! Where the total node counter ends nodes
integer, intent(in) :: end_n_in(3)

integer :: fn=1        ! CGNS file index number
integer :: ier         ! CGNS error status
integer :: base=1      ! base number
integer :: zone=1      ! zone number
integer :: nnodes      ! Number of nodes in this processor
integer :: sol =1      ! solution number
integer :: field       ! section number
integer(cgsize_t) :: sizes(3,3)  ! Sizes

! Convert input to right data type
integer(cgsize_t) :: start_n(3)  ! Where the total node counter starts nodes
integer(cgsize_t) :: end_n(3)  ! Where the total node counter ends nodes

! Building the lcoal mesh
integer :: i,j,k
real(rprec), dimension(nx,ny,nz) :: xyz

!  ! Set the parallel communicator
!  call cgp_mpi_comm_f(cgnsParallelComm, ierr)

! Convert types such that CGNS libraries can handle the input
start_n(1) = int(start_n_in(1), cgsize_t)
start_n(2) = int(start_n_in(2), cgsize_t)
start_n(3) = int(start_n_in(3), cgsize_t)
end_n(1) = int(end_n_in(1), cgsize_t)
end_n(2) = int(end_n_in(2), cgsize_t)
end_n(3) = int(end_n_in(3), cgsize_t)

! The total number of nodes in this processor
nnodes = nx*ny*nz

! Sizes, used to create zone
sizes(:,1) = (/int(nx, cgsize_t),int(ny, cgsize_t),int(nz_tot, cgsize_t)/)
sizes(:,2) = (/int(nx-1, cgsize_t),int(ny-1, cgsize_t),int(nz_tot-1, cgsize_t)/)
sizes(:,3) = (/int(0, cgsize_t) , int(0, cgsize_t), int(0, cgsize_t)/)

! Open CGNS file
call cgp_open_f(file_name, CG_MODE_WRITE, fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write base
call cg_base_write_f(fn, 'Base', 3, 3, base, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write zone
call cg_zone_write_f(fn, base, 'Zone', sizes, Structured, zone, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write print info to screen
if (coord .eq. 0) then
    write(*,*) 'Writing, ', file_name
end if

! Create data nodes for coordinates
call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateX', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateY', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

call cgp_coord_write_f(fn, base, zone, RealDouble, 'CoordinateZ', nnodes, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! This is done for the 3 dimensions x,y and z
! It writes the coordinates
! Create grid points
do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = xin(i)
end do
end do
end do
write(*,*) "HERE 0.8"

call cgp_coord_write_data_f(fn, base, zone, 1, start_n, end_n, %VAL(0), ier)
write(*,*) "HERE 0.85"
if (ier .ne. CG_OK) call cgp_error_exit_f
write(*,*) "HERE 0.9"

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = yin(j)
end do
end do
end do
call cgp_coord_write_data_f(fn, base, zone, 2, start_n, end_n, %VAL(0), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f
write(*,*) "HERE 1.0"

! Write out the queued coordinate data
!  call cgp_queue_flush_f(ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f
!  call cgp_queue_set_f(0, ier)

! Write the coordinate data in parallel to the queue
!  call cgp_queue_set_f(1, ier)
!  if (ier .ne. CG_OK) call cgp_error_exit_f

do k = 1, nz
do j = 1, ny
do i = 1, nx
    xyz(i,j,k) = zin(k)
end do
end do
end do
write(*,*) "HERE 1.1"

call cgp_coord_write_data_f(fn, base, zone, 3, start_n, end_n, %VAL(0), ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Create a centered solution
call cg_sol_write_f(fn, base, zone, 'Solution', Vertex, sol, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

! Write the solution
do i = 1, num_fields
    call cgp_field_write_f(fn, base, zone, sol, RealDouble, fieldNames(i),     &
                           field, ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

    call cgp_field_write_data_f(fn, base, zone, sol, field, start_n, end_n,    &
                                %VAL(0), ier)
    if (ier .ne. CG_OK) call cgp_error_exit_f

end do

! Close the file
call cgp_close_f(fn, ier)
if (ier .ne. CG_OK) call cgp_error_exit_f

write(*,*) "end of write_null_cgns"

end subroutine write_null_cgns
#endif
#endif

!*******************************************************************************
subroutine output_loop()
!*******************************************************************************
!
!  This subroutine is called every time step and acts as a driver for
!  computing statistics and outputing instantaneous data. No actual
!  calculations are performed here.
!
use param, only : jt_total, dt
use param, only : checkpoint_data, checkpoint_nskip
use param, only : tavg_calc, tavg_nstart, tavg_nend, tavg_nskip
use param, only : point_calc, point_nstart, point_nend, point_nskip
use param, only : domain_calc, domain_nstart, domain_nend, domain_nskip
use param, only : xplane_calc, xplane_nstart, xplane_nend, xplane_nskip
use param, only : yplane_calc, yplane_nstart, yplane_nend, yplane_nskip
use param, only : zplane_calc, zplane_nstart, zplane_nend, zplane_nskip
implicit none

! Determine if we are to checkpoint intermediate times
if( checkpoint_data ) then
    ! Now check if data should be checkpointed this time step
    if ( modulo (jt_total, checkpoint_nskip) == 0) call checkpoint()
end if

!  Determine if time summations are to be calculated
if (tavg_calc) then
    ! Are we between the start and stop timesteps?
    if ((jt_total >= tavg_nstart).and.(jt_total <= tavg_nend)) then
        ! Every timestep (between nstart and nend), add to tavg%dt
        tavg%dt = tavg%dt + dt

        ! Are we at the beginning or a multiple of nstart?
        if ( mod(jt_total-tavg_nstart,tavg_nskip)==0 ) then
            ! Check if we have initialized tavg
            if (.not.tavg%initialized) then
                if (coord == 0) then
                    write(*,*) '-------------------------------'
                    write(*,"(1a,i9,1a,i9)")                                   &
                        'Starting running time summation from ',               &
                        tavg_nstart, ' to ', tavg_nend
                    write(*,*) '-------------------------------'
                end if

                call tavg%init()
            else
                call tavg%compute()
            end if
        end if
    end if
end if

!  Determine if instantaneous point velocities are to be recorded
if(point_calc) then
    if (jt_total >= point_nstart .and. jt_total <= point_nend .and.            &
        ( mod(jt_total-point_nstart,point_nskip)==0) ) then
        if (jt_total == point_nstart) then
            if (coord == 0) then
                write(*,*) '-------------------------------'
                write(*,"(1a,i9,1a,i9)")                                       &
                    'Writing instantaneous point velocities from ',            &
                    point_nstart, ' to ', point_nend
                write(*,"(1a,i9)") 'Iteration skip:', point_nskip
                write(*,*) '-------------------------------'
            end if
        end if
        call inst_write(1)
    end if
end if

!  Determine if instantaneous domain velocities are to be recorded
if(domain_calc) then
    if (jt_total >= domain_nstart .and. jt_total <= domain_nend .and.          &
        ( mod(jt_total-domain_nstart,domain_nskip)==0) ) then
        if (jt_total == domain_nstart) then
            if (coord == 0) then
                write(*,*) '-------------------------------'
                write(*,"(1a,i9,1a,i9)")                                       &
                    'Writing instantaneous domain velocities from ',           &
                    domain_nstart, ' to ', domain_nend
                write(*,"(1a,i9)") 'Iteration skip:', domain_nskip
                write(*,*) '-------------------------------'
            end if

        end if
        call inst_write(2)
    end if
end if

!  Determine if instantaneous x-plane velocities are to be recorded
if(xplane_calc) then
    if (jt_total >= xplane_nstart .and. jt_total <= xplane_nend .and.          &
        ( mod(jt_total-xplane_nstart,xplane_nskip)==0) ) then
    if (jt_total == xplane_nstart) then
        if (coord == 0) then
            write(*,*) '-------------------------------'
            write(*,"(1a,i9,1a,i9)")                                           &
                'Writing instantaneous x-plane velocities from ',              &
                xplane_nstart, ' to ', xplane_nend
            write(*,"(1a,i9)") 'Iteration skip:', xplane_nskip
            write(*,*) '-------------------------------'
            end if
        end if

        call inst_write(3)
    end if
end if

!  Determine if instantaneous y-plane velocities are to be recorded
if(yplane_calc) then
    if (jt_total >= yplane_nstart .and. jt_total <= yplane_nend .and.          &
        ( mod(jt_total-yplane_nstart,yplane_nskip)==0) ) then
        if (jt_total == yplane_nstart) then
            if (coord == 0) then
                write(*,*) '-------------------------------'
                write(*,"(1a,i9,1a,i9)")                                       &
                    'Writing instantaneous y-plane velocities from ',          &
                    yplane_nstart, ' to ', yplane_nend
                write(*,"(1a,i9)") 'Iteration skip:', yplane_nskip
                write(*,*) '-------------------------------'
            end if
        end if

        call inst_write(4)
    end if
end if

!  Determine if instantaneous z-plane velocities are to be recorded
if(zplane_calc) then
    if (jt_total >= zplane_nstart .and. jt_total <= zplane_nend .and.          &
        ( mod(jt_total-zplane_nstart,zplane_nskip)==0) ) then
        if (jt_total == zplane_nstart) then
            if (coord == 0) then
                write(*,*) '-------------------------------'
                write(*,"(1a,i9,1a,i9)")                                       &
                    'Writing instantaneous z-plane velocities from ',          &
                    zplane_nstart, ' to ', zplane_nend
                write(*,"(1a,i9)") 'Iteration skip:', zplane_nskip
                write(*,*) '-------------------------------'
            end if
        end if

        call inst_write(5)
    end if
end if

end subroutine output_loop

!*******************************************************************************
subroutine inst_write(itype)
!*******************************************************************************
!
! This subroutine is used to write all of the instantaneous data from
! lesgo to file. The types of data written are:
!
!   points   : itype=1
!   domain   : itype=2
!   x-planes : itype=3
!   y-planes : itype=4
!   z-planes : itype=5
!
! For the points and planar data, this subroutine writes using the
! locations specfied from the param module.
! If additional instantenous values are
! desired to be written, they should be done so using this subroutine.
!
use functions, only : linear_interp, trilinear_interp, interp_to_uv_grid
use param, only : point_nloc, point_loc
use param, only : xplane_nloc, xplane_loc
use param, only : yplane_nloc, yplane_loc
use param, only : zplane_nloc, zplane_loc
use param, only : dx, dy, dz
use param, only : write_endian
use grid_m, only : grid
use sim_param, only : u, v, w, p
use sim_param, only : dwdy, dwdx, dvdx, dudy
use functions, only : interp_to_w_grid

use stat_defs, only : point, xplane, yplane
#ifdef PPMPI
use stat_defs, only : zplane
use param, only : ny, nz, nproc
#endif
#ifdef PPLVLSET
use level_set_base, only : phi
use sim_param, only : fx, fy, fz, fxa, fya, fza
#endif
#ifdef PPSCALARS
use scalars, only : theta
#endif

implicit none

integer, intent(in) :: itype
character (64) :: fname
integer :: n, i, j, k
real(rprec), allocatable, dimension(:,:,:) :: ui, vi, wi,w_uv
#if defined(PPLES_GPU)
integer :: p_i0, p_j0, p_k0, p_i1, p_j1, p_k1
integer, allocatable, dimension(:) :: point_active, point_i, point_j, point_k
real(rprec), allocatable, dimension(:) :: point_xd, point_yd, point_zd
real(rprec), allocatable, dimension(:,:) :: point_vals
integer :: z_k0, z_k1
integer, allocatable, dimension(:) :: zplane_active, zplane_k
real(rprec), allocatable, dimension(:) :: zplane_l
real(rprec), allocatable, dimension(:,:,:) :: z_ui, z_vi, z_wi
integer :: x_i0, x_i1, y_j0, y_j1
integer, allocatable, dimension(:) :: xplane_i, yplane_j
real(rprec), allocatable, dimension(:) :: xplane_l, yplane_l
real(rprec), allocatable, dimension(:,:,:) :: x_ui, x_vi, x_wi
real(rprec), allocatable, dimension(:,:,:) :: y_ui, y_vi, y_wi
real(rprec) :: px, py, pz
real(rprec) :: u1, u2, u3, u4, u5, u6
real(rprec) :: w1, w2, w3, w4, wuv0, wuv1, zl
#endif
real(rprec), pointer, dimension(:) :: x, y, z, zw
! Vorticity
real(rprec), dimension (:,:,:), allocatable :: vortx, vorty, vortz

! Pressure
real(rprec), dimension(:,:,:), allocatable :: pres_real
#ifndef PPCGNS
character(64) :: bin_ext

#ifdef PPLVLSET
real(rprec), allocatable, dimension(:,:,:) :: fx_tot, fy_tot, fz_tot
#endif

#ifdef PPMPI
call string_splice(bin_ext, '.c', coord, '.bin')
#else
bin_ext = '.bin'
#endif
#endif

! Nullify pointers
nullify(x,y,z,zw)

! Set grid pointers
x => grid % x
y => grid % y
z => grid % z
zw => grid % zw

#if defined(PPLES_GPU)
! Fast path for single-rank point output: sample the resident device fields
! into a small buffer and transfer only O(point_nloc) values. Multi-rank point
! output still uses the host path below until a device halo-aware w_uv sampler
! is validated.
#ifdef PPMPI
if (itype == 1 .and. nproc == 1) then
#else
if (itype == 1) then
#endif
    allocate(point_active(point_nloc), point_i(point_nloc),                   &
             point_j(point_nloc), point_k(point_nloc))
    allocate(point_xd(point_nloc), point_yd(point_nloc),                     &
             point_zd(point_nloc))
    allocate(point_vals(3, point_nloc))

    point_active = 0
    point_i = 1
    point_j = 1
    point_k = 1
    point_xd = 0._rprec
    point_yd = 0._rprec
    point_zd = 0._rprec
    point_vals = 0._rprec

    do n = 1, point_nloc
#ifdef PPMPI
        if (point(n) % coord == 0) point_active(n) = 1
#else
        point_active(n) = 1
#endif
        point_i(n) = point(n) % istart
        point_j(n) = point(n) % jstart
        point_k(n) = point(n) % kstart
        point_xd(n) = point(n) % xdiff
        point_yd(n) = point(n) % ydiff
        point_zd(n) = point(n) % zdiff
    end do

    !$acc wait(1)
    !$acc data copyin(point_active, point_i, point_j, point_k,                &
    !$acc&            point_xd, point_yd, point_zd) copyout(point_vals)
    !$acc parallel loop default(present) private(p_i0, p_j0, p_k0, p_i1,     &
    !$acc&            p_j1, p_k1, px, py, pz, u1, u2, u3, u4, u5, u6,        &
    !$acc&            w1, w2, w3, w4)
    do n = 1, point_nloc
        if (point_active(n) /= 0) then
            p_i0 = point_i(n)
            p_j0 = point_j(n)
            p_k0 = point_k(n)
            p_i1 = p_i0 + 1
            p_j1 = p_j0 + 1
            p_k1 = p_k0 + 1
            if (p_i1 > nx) p_i1 = 1
            if (p_j1 > ny) p_j1 = 1

            px = point_xd(n)
            py = point_yd(n)
            pz = point_zd(n)

            u1 = u(p_i0,p_j0,p_k0) + px * (u(p_i1,p_j0,p_k0)                 &
                - u(p_i0,p_j0,p_k0)) / dx
            u2 = u(p_i0,p_j1,p_k0) + px * (u(p_i1,p_j1,p_k0)                 &
                - u(p_i0,p_j1,p_k0)) / dx
            u3 = u(p_i0,p_j0,p_k1) + px * (u(p_i1,p_j0,p_k1)                 &
                - u(p_i0,p_j0,p_k1)) / dx
            u4 = u(p_i0,p_j1,p_k1) + px * (u(p_i1,p_j1,p_k1)                 &
                - u(p_i0,p_j1,p_k1)) / dx
            u5 = u1 + py * (u2 - u1) / dy
            u6 = u3 + py * (u4 - u3) / dy
            point_vals(1,n) = u5 + pz * (u6 - u5) / dz

            u1 = v(p_i0,p_j0,p_k0) + px * (v(p_i1,p_j0,p_k0)                 &
                - v(p_i0,p_j0,p_k0)) / dx
            u2 = v(p_i0,p_j1,p_k0) + px * (v(p_i1,p_j1,p_k0)                 &
                - v(p_i0,p_j1,p_k0)) / dx
            u3 = v(p_i0,p_j0,p_k1) + px * (v(p_i1,p_j0,p_k1)                 &
                - v(p_i0,p_j0,p_k1)) / dx
            u4 = v(p_i0,p_j1,p_k1) + px * (v(p_i1,p_j1,p_k1)                 &
                - v(p_i0,p_j1,p_k1)) / dx
            u5 = u1 + py * (u2 - u1) / dy
            u6 = u3 + py * (u4 - u3) / dy
            point_vals(2,n) = u5 + pz * (u6 - u5) / dz

            w1 = 0.5_rprec * (w(p_i0,p_j0,p_k0+1) + w(p_i0,p_j0,p_k0))
            w2 = 0.5_rprec * (w(p_i1,p_j0,p_k0+1) + w(p_i1,p_j0,p_k0))
            w3 = 0.5_rprec * (w(p_i0,p_j1,p_k0+1) + w(p_i0,p_j1,p_k0))
            w4 = 0.5_rprec * (w(p_i1,p_j1,p_k0+1) + w(p_i1,p_j1,p_k0))
            u1 = w1 + px * (w2 - w1) / dx
            u2 = w3 + px * (w4 - w3) / dx

            if (p_k1 < nz) then
                w1 = 0.5_rprec * (w(p_i0,p_j0,p_k1+1) + w(p_i0,p_j0,p_k1))
                w2 = 0.5_rprec * (w(p_i1,p_j0,p_k1+1) + w(p_i1,p_j0,p_k1))
                w3 = 0.5_rprec * (w(p_i0,p_j1,p_k1+1) + w(p_i0,p_j1,p_k1))
                w4 = 0.5_rprec * (w(p_i1,p_j1,p_k1+1) + w(p_i1,p_j1,p_k1))
                u3 = w1 + px * (w2 - w1) / dx
                u4 = w3 + px * (w4 - w3) / dx
                u6 = u3 + py * (u4 - u3) / dy
            else
                u6 = u1 + py * (u2 - u1) / dy
            end if
            u5 = u1 + py * (u2 - u1) / dy
            point_vals(3,n) = u5 + pz * (u6 - u5) / dz
        end if
    end do
    !$acc end data

    do n = 1, point_nloc
        call string_splice(fname, path // 'output/vel.x-', point_loc(n)%xyz(1),&
            '.y-', point_loc(n)%xyz(2), '.z-', point_loc(n)%xyz(3), '.dat')
#ifdef PPMPI
        if(point(n) % coord == 0) then
#endif
            open(unit=13, position="append", file=fname)
            write(13,*) total_time, point_vals(1,n), point_vals(2,n),          &
                point_vals(3,n)
            close(13)
#ifdef PPMPI
        end if
#endif
    end do

    deallocate(point_active, point_i, point_j, point_k, point_xd, point_yd,    &
               point_zd, point_vals)
    nullify(x,y,z,zw)
    return
end if

#if defined(PPMPI) && !defined(PPCGNS)
! Fast path for single-rank x/y-plane velocity output. This mirrors the host
! interpolation formulas while transferring only plane-sized buffers.
if ((itype == 3 .or. itype == 4) .and. nproc == 1) then
    if (itype == 3) then
        allocate(xplane_i(xplane_nloc), xplane_l(xplane_nloc))
        allocate(x_ui(ny,nz,xplane_nloc), x_vi(ny,nz,xplane_nloc),            &
                 x_wi(ny,nz,xplane_nloc))
        x_ui = 0._rprec
        x_vi = 0._rprec
        x_wi = 0._rprec
        do n = 1, xplane_nloc
            xplane_i(n) = xplane(n) % istart
            xplane_l(n) = xplane(n) % ldiff
        end do

        !$acc wait(1)
        !$acc data copyin(xplane_i, xplane_l) copyout(x_ui, x_vi, x_wi)
        !$acc parallel loop collapse(3) default(present) private(x_i0, x_i1,  &
        !$acc&            zl, wuv0, wuv1)
        do n = 1, xplane_nloc
        do k = 1, nz
        do j = 1, ny
            x_i0 = xplane_i(n)
            x_i1 = x_i0 + 1
            zl = xplane_l(n)
            x_ui(j,k,n) = u(x_i0,j,k) + zl * (u(x_i1,j,k) - u(x_i0,j,k)) / dx
            x_vi(j,k,n) = v(x_i0,j,k) + zl * (v(x_i1,j,k) - v(x_i0,j,k)) / dx
            if (k < nz) then
                wuv0 = 0.5_rprec * (w(x_i0,j,k+1) + w(x_i0,j,k))
                wuv1 = 0.5_rprec * (w(x_i1,j,k+1) + w(x_i1,j,k))
            else
                wuv0 = 0.5_rprec * (w(x_i0,j,k) + w(x_i0,j,k-1))
                wuv1 = 0.5_rprec * (w(x_i1,j,k) + w(x_i1,j,k-1))
            end if
            x_wi(j,k,n) = wuv0 + zl * (wuv1 - wuv0) / dx
        end do
        end do
        end do
        !$acc end data

        do n = 1, xplane_nloc
            call string_splice(fname, path // 'output/vel.x-',                &
                 xplane_loc(n), '.', jt_total)
            call string_concat(fname, bin_ext)
            open(unit=13,file=fname,form='unformatted',convert=write_endian,  &
                access='direct',recl=ny*nz*rprec)
            write(13,rec=1) x_ui(1:ny,1:nz,n)
            write(13,rec=2) x_vi(1:ny,1:nz,n)
            write(13,rec=3) x_wi(1:ny,1:nz,n)
            close(13)
        end do

        deallocate(xplane_i, xplane_l, x_ui, x_vi, x_wi)
    else
        allocate(yplane_j(yplane_nloc), yplane_l(yplane_nloc))
        allocate(y_ui(nx,nz,yplane_nloc), y_vi(nx,nz,yplane_nloc),            &
                 y_wi(nx,nz,yplane_nloc))
        y_ui = 0._rprec
        y_vi = 0._rprec
        y_wi = 0._rprec
        do n = 1, yplane_nloc
            yplane_j(n) = yplane(n) % istart
            yplane_l(n) = yplane(n) % ldiff
        end do

        !$acc wait(1)
        !$acc data copyin(yplane_j, yplane_l) copyout(y_ui, y_vi, y_wi)
        !$acc parallel loop collapse(3) default(present) private(y_j0, y_j1,  &
        !$acc&            zl, wuv0, wuv1)
        do n = 1, yplane_nloc
        do k = 1, nz
        do i = 1, nx
            y_j0 = yplane_j(n)
            y_j1 = y_j0 + 1
            zl = yplane_l(n)
            y_ui(i,k,n) = u(i,y_j0,k) + zl * (u(i,y_j1,k) - u(i,y_j0,k)) / dy
            y_vi(i,k,n) = v(i,y_j0,k) + zl * (v(i,y_j1,k) - v(i,y_j0,k)) / dy
            if (k < nz) then
                wuv0 = 0.5_rprec * (w(i,y_j0,k+1) + w(i,y_j0,k))
                wuv1 = 0.5_rprec * (w(i,y_j1,k+1) + w(i,y_j1,k))
            else
                wuv0 = 0.5_rprec * (w(i,y_j0,k) + w(i,y_j0,k-1))
                wuv1 = 0.5_rprec * (w(i,y_j1,k) + w(i,y_j1,k-1))
            end if
            y_wi(i,k,n) = wuv0 + zl * (wuv1 - wuv0) / dy
        end do
        end do
        end do
        !$acc end data

        do n = 1, yplane_nloc
            call string_splice(fname, path // 'output/vel.y-',                &
                 yplane_loc(n), '.', jt_total)
            call string_concat(fname, bin_ext)
            open(unit=13,file=fname,form='unformatted',convert=write_endian,  &
                access='direct',recl=nx*nz*rprec)
            write(13,rec=1) y_ui(1:nx,1:nz,n)
            write(13,rec=2) y_vi(1:nx,1:nz,n)
            write(13,rec=3) y_wi(1:nx,1:nz,n)
            close(13)
        end do

        deallocate(yplane_j, yplane_l, y_ui, y_vi, y_wi)
    end if

    nullify(x,y,z,zw)
    return
end if

! Fast path for single-rank z-plane velocity output. It copies only the
! requested plane buffers to host instead of refreshing full u/v/w/p/gradient
! fields. Multi-rank needs a halo-aware device sampler before using this path.
if (itype == 5 .and. nproc == 1) then
    allocate(zplane_active(zplane_nloc), zplane_k(zplane_nloc),               &
             zplane_l(zplane_nloc))
    allocate(z_ui(nx,ny,zplane_nloc), z_vi(nx,ny,zplane_nloc),                &
             z_wi(nx,ny,zplane_nloc))

    zplane_active = 0
    zplane_k = 1
    zplane_l = 0._rprec
    z_ui = 0._rprec
    z_vi = 0._rprec
    z_wi = 0._rprec

    do k = 1, zplane_nloc
        if (zplane(k) % coord == 0) zplane_active(k) = 1
        zplane_k(k) = zplane(k) % istart
        zplane_l(k) = zplane(k) % ldiff
    end do

    !$acc wait(1)
    !$acc data copyin(zplane_active, zplane_k, zplane_l) copyout(z_ui, z_vi, z_wi)
    !$acc parallel loop collapse(3) default(present) private(z_k0, z_k1,      &
    !$acc&            zl, wuv0, wuv1)
    do k = 1, zplane_nloc
    do j = 1, ny
    do i = 1, nx
        if (zplane_active(k) /= 0) then
            z_k0 = zplane_k(k)
            z_k1 = z_k0 + 1
            zl = zplane_l(k)
            z_ui(i,j,k) = u(i,j,z_k0) + zl * (u(i,j,z_k1) - u(i,j,z_k0)) / dz
            z_vi(i,j,k) = v(i,j,z_k0) + zl * (v(i,j,z_k1) - v(i,j,z_k0)) / dz

            wuv0 = 0.5_rprec * (w(i,j,z_k0+1) + w(i,j,z_k0))
            if (z_k1 < nz) then
                wuv1 = 0.5_rprec * (w(i,j,z_k1+1) + w(i,j,z_k1))
            else
                wuv1 = wuv0
            end if
            z_wi(i,j,k) = wuv0 + zl * (wuv1 - wuv0) / dz
        end if
    end do
    end do
    end do
    !$acc end data

    do k = 1, zplane_nloc
        call string_splice(fname, path // 'output/vel.z-',                    &
                zplane_loc(k), '.', jt_total)
        if(zplane(k) % coord == 0) then
            call string_concat(fname, bin_ext)
            open(unit=13,file=fname,form='unformatted',convert=write_endian,  &
                            access='direct',recl=nx*ny*1*rprec)
            write(13,rec=1) z_ui(1:nx,1:ny,k)
            write(13,rec=2) z_vi(1:nx,1:ny,k)
            write(13,rec=3) z_wi(1:nx,1:ny,k)
            close(13)
        end if
    end do

    deallocate(zplane_active, zplane_k, zplane_l, z_ui, z_vi, z_wi)
    nullify(x,y,z,zw)
    return
end if
#endif

! Other instantaneous output modes write domain or plane-sized host files, so
! keep the existing full refresh before their host-side interpolation/writes.
!$acc wait(1)
!$acc update self(u, v, w, p, dwdy, dwdx, dvdx, dudy)
#endif

!  Allocate space for the interpolated w values
allocate(w_uv(nx,ny,lbz:nz))

!  Make sure w has been interpolated to uv-grid
w_uv = interp_to_uv_grid(w(1:nx,1:ny,lbz:nz), lbz)

!  Instantaneous velocity sampled at point
if(itype==1) then
    do n = 1, point_nloc
        ! Common file name for all output types
        call string_splice(fname, path // 'output/vel.x-', point_loc(n)%xyz(1),&
            '.y-', point_loc(n)%xyz(2), '.z-', point_loc(n)%xyz(3), '.dat')

#ifdef PPMPI
        if(point(n) % coord == coord) then
#endif
            open(unit=13, position="append", file=fname)
            write(13,*) total_time,                                            &
            trilinear_interp(u(1:nx,1:ny,lbz:nz), lbz, point_loc(n)%xyz),      &
            trilinear_interp(v(1:nx,1:ny,lbz:nz), lbz, point_loc(n)%xyz),      &
            trilinear_interp(w_uv(1:nx,1:ny,lbz:nz), lbz, point_loc(n)%xyz)
            close(13)
#ifdef PPMPI
        end if
#endif
    end do

!  Instantaneous write for entire domain
elseif(itype==2) then
    ! Common file name for all output types
    call string_splice(fname, path //'output/vel.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                 &
        3, (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),                        &
        (/ u(1:nx,1:ny,1:(nz-nz_end)), v(1:nx,1:ny,1:(nz-nz_end)),             &
         w_uv(1:nx,1:ny,1:(nz-nz_end)) /) )
#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=nx*ny*nz*rprec)
    write(13,rec=1) u(:nx,:ny,1:nz)
    write(13,rec=2) v(:nx,:ny,1:nz)
    write(13,rec=3) w_uv(:nx,:ny,1:nz)
    close(13)
#endif

    ! Compute vorticity
    allocate(vortx(nx,ny,lbz:nz), vorty(nx,ny,lbz:nz), vortz(nx,ny,lbz:nz))
    vortx(1:nx,1:ny,lbz:nz) = 0._rprec
    vorty(1:nx,1:ny,lbz:nz) = 0._rprec
    vortz(1:nx,1:ny,lbz:nz) = 0._rprec

    ! Use vorticityx as an intermediate step for performing uv-w interpolation
    ! Vorticity is written in w grid
    vortx(1:nx,1:ny,lbz:nz) = dvdx(1:nx,1:ny,lbz:nz) - dudy(1:nx,1:ny,lbz:nz)
    vortz(1:nx,1:ny,lbz:nz) = interp_to_w_grid( vortx(1:nx,1:ny,lbz:nz), lbz)
    vortx(1:nx,1:ny,lbz:nz) = dwdy(1:nx,1:ny,lbz:nz) - dvdz(1:nx,1:ny,lbz:nz)
    vorty(1:nx,1:ny,lbz:nz) = dudz(1:nx,1:ny,lbz:nz) - dwdx(1:nx,1:ny,lbz:nz)

    if (coord == 0) then
        vortz(1:nx,1:ny, 1) = 0._rprec
    end if

    ! Common file name for all output types
    call string_splice(fname, path //'output/vort.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname,nx,ny, nz - nz_end, nz_tot,                 &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nx) , y(1:ny) , zw(1:(nz-nz_end) ),                                &
        3, (/ 'VorticityX', 'VorticityY', 'VorticityZ' /),                     &
        (/ vortx(1:nx,1:ny,1:(nz-nz_end)), vorty(1:nx,1:ny,1:(nz-nz_end)),     &
        vortz(1:nx,1:ny,1:(nz-nz_end)) /) )

#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=nx*ny*nz*rprec)
    write(13,rec=1) vortx(:nx,:ny,1:nz)
    write(13,rec=2) vorty(:nx,:ny,1:nz)
    write(13,rec=3) vortz(:nx,:ny,1:nz)
    close(13)
#endif

    deallocate(vortx, vorty, vortz)

    ! Compute pressure
    allocate(pres_real(nx,ny,lbz:nz))
    pres_real(1:nx,1:ny,lbz:nz) = 0._rprec

    ! Calculate real pressure
    pres_real(1:nx,1:ny,lbz:nz) = p(1:nx,1:ny,lbz:nz)                          &
        - 0.5 * ( u(1:nx,1:ny,lbz:nz)**2                                       &
        + interp_to_uv_grid( w(1:nx,1:ny,lbz:nz), lbz)**2                      &
        + v(1:nx,1:ny,lbz:nz)**2 )

    ! Common file name for all output types
    call string_splice(fname, path //'output/pres.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
        (/ 1, 1,   (nz-1)*coord + 1 /),                                        &
        (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                           &
        x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                 &
        1, (/ 'Pressure' /), (/ pres_real(1:nx,1:ny,1:(nz-nz_end)) /) )

#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
        access='direct', recl=nx*ny*nz*rprec)
    write(13,rec=1) pres_real(:nx,:ny,1:nz)
    close(13)
#endif

     deallocate(pres_real)

#ifdef PPSCALARS
    ! Common file name for all output types
    call string_splice(fname, path //'output/theta.', jt_total)
#if defined(PPCGNS) && defined(PPMPI)
    ! Write CGNS Output
    call string_concat(fname, '.cgns')
    call write_parallel_cgns(fname, nx, ny, nz - nz_end, nz_tot,               &
     (/ 1, 1,   (nz-1)*coord + 1 /),                                           &
     (/ nx, ny, (nz-1)*(coord+1) + 1 - nz_end /),                              &
     x(1:nx) , y(1:ny) , z(1:(nz-nz_end) ),                                    &
     1, (/ 'Theta' /), (/ theta(1:nx,1:ny,1:(nz-nz_end)) /) )
#else
    ! Write binary Output
    call string_concat(fname, bin_ext)
    open(unit=13, file=fname, form='unformatted', convert=write_endian,        &
     access='direct', recl=nx*ny*nz*rprec)
    write(13,rec=1) theta(:nx,:ny,1:nz)
    close(13)
#endif
#endif

!  Write instantaneous x-plane values
elseif(itype==3) then

    allocate(ui(1,ny,nz), vi(1,ny,nz), wi(1,ny,nz))

    !  Loop over all xplane locations
    do i = 1, xplane_nloc
        do k = 1, nz
            do j = 1, ny
                ui(1,j,k) = linear_interp(u(xplane(i) % istart,j,k),    &
                     u(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)
                vi(1,j,k) = linear_interp(v(xplane(i) % istart,j,k),    &
                     v(xplane(i) % istart+1,j,k), dx, xplane(i) % ldiff)
                wi(1,j,k) = linear_interp(w_uv(xplane(i) % istart,j,k), &
                     w_uv(xplane(i) % istart+1,j,k), dx, &
                     xplane(i) % ldiff)
            end do
        end do

        ! Common file name portion for all output types
        call string_splice(fname, path // 'output/vel.x-', xplane_loc(i), '.', jt_total)

#if defined(PPCGNS) && defined(PPMPI)
        ! Write CGNS Output
        call string_concat(fname, '.cgns')
        call write_parallel_cgns (fname,1,ny, nz - nz_end, nz_tot,     &
                        (/ 1, 1,   (nz-1)*coord + 1 /),                &
                        (/ 1, ny, (nz-1)*(coord+1) + 1 - nz_end /),    &
                    xplane_loc(i:i) , y(1:ny) , z(1:(nz-nz_end) ),     &
              3, (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),          &
              (/ ui(1,1:ny,1:(nz-nz_end)), vi(1,1:ny,1:(nz-nz_end)),   &
                 wi(1,1:ny,1:(nz-nz_end)) /) )

#else
        ! Write binary output
        call string_concat(fname, bin_ext)
        open(unit=13,file=fname,form='unformatted',convert=write_endian, access='direct',recl=ny*nz*rprec)
        write(13,rec=1) ui
        write(13,rec=2) vi
        write(13,rec=3) wi
        close(13)
#endif
    end do

    deallocate(ui,vi,wi)

!  Write instantaneous y-plane values
elseif(itype==4) then

    allocate(ui(nx,1,nz), vi(nx,1,nz), wi(nx,1,nz))

    !  Loop over all yplane locations
    do j = 1, yplane_nloc
        do k = 1, nz
            do i = 1, nx

                ui(i,1,k) = linear_interp(u(i,yplane(j) % istart,k),           &
                     u(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
                vi(i,1,k) = linear_interp(v(i,yplane(j) % istart,k),           &
                     v(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
                wi(i,1,k) = linear_interp(w_uv(i,yplane(j) % istart,k),        &
                     w_uv(i,yplane(j) % istart+1,k), dy, yplane(j) % ldiff)
            end do
        end do

        ! Common file name portion for all output types
        call string_splice(fname, path // 'output/vel.y-', yplane_loc(j), '.', &
             jt_total)

#if defined(PPCGNS) && defined(PPMPI)
        call string_concat(fname, '.cgns')
        call write_parallel_cgns (fname,nx,1, nz - nz_end, nz_tot,             &
            (/ 1, 1,   (nz-1)*coord + 1 /),                                    &
            (/ nx, 1, (nz-1)*(coord+1) + 1 - nz_end /),                        &
            x(1:nx) , yplane_loc(j:j) , z(1:(nz-nz_end) ),                     &
            3, (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),                    &
            (/ ui(1:nx,1,1:(nz-nz_end)), vi(1:nx,1,1:(nz-nz_end)),             &
            wi(1:nx,1,1:(nz-nz_end)) /) )
#else
        ! Write binary output
        call string_concat(fname, bin_ext)
        open(unit=13,file=fname,form='unformatted',convert=write_endian, access='direct',recl=nx*nz*rprec)
        write(13,rec=1) ui
        write(13,rec=2) vi
        write(13,rec=3) wi
        close(13)
#endif

    end do

    deallocate(ui,vi,wi)

!  Write instantaneous z-plane values
elseif (itype==5) then

    allocate(ui(nx,ny,1), vi(nx,ny,1), wi(nx,ny,1))

    !  Loop over all zplane locations
    do k = 1, zplane_nloc
        ! Common file name portion for all output types
        call string_splice(fname, path // 'output/vel.z-',                     &
                zplane_loc(k), '.', jt_total)

#ifdef PPCGNS
        call string_concat(fname, '.cgns')
#endif

#ifdef PPMPI
        if(zplane(k) % coord == coord) then
            do j = 1, Ny
                do i = 1, Nx
                    ui(i,j,1) = linear_interp(u(i,j,zplane(k) % istart),       &
                         u(i,j,zplane(k) % istart+1), dz, zplane(k) % ldiff)
                    vi(i,j,1) = linear_interp(v(i,j,zplane(k) % istart),       &
                         v(i,j,zplane(k) % istart+1), dz, zplane(k) % ldiff)
                    wi(i,j,1) = linear_interp(w_uv(i,j,zplane(k) % istart),    &
                         w_uv(i,j,zplane(k) % istart+1), dz, zplane(k) % ldiff)
                end do
            end do

#ifdef PPCGNS
            call warn("inst_write","Z plane writting is currently disabled.")
!            ! Write CGNS Data
!            ! Only the processor with data writes, the other one is written
!            ! using null arguments with 'write_null_cgns'
!            call write_parallel_cgns (fname ,nx, ny, 1, 1,                     &
!                (/ 1, 1,   1 /),                                               &
!                (/ nx, ny, 1 /),                                               &
!                x(1:nx) , y(1:ny) , zplane_loc(k:k), 3,                        &
!                (/ 'VelocityX', 'VelocityY', 'VelocityZ' /),                   &
!                (/ ui(1:nx,1:ny,1), vi(1:nx,1:ny,1), wi(1:nx,1:ny,1) /) )
#else
            call string_concat(fname, bin_ext)
            open(unit=13,file=fname,form='unformatted',convert=write_endian,   &
                            access='direct',recl=nx*ny*1*rprec)
            write(13,rec=1) ui(1:nx,1:ny,1)
            write(13,rec=2) vi(1:nx,1:ny,1)
            write(13,rec=3) wi(1:nx,1:ny,1)
            close(13)
#endif
!
! #ifdef PPMPI
!         else
! #ifdef PPCGNS
!            write(*,*) "At write_null_cgns"
!            call write_null_cgns (fname ,nx, ny, 1, 1,                         &
!            (/ 1, 1,   1 /),                                                   &
!            (/ nx, ny, 1 /),                                                   &
!            x(1:nx) , y(1:ny) , zplane_loc(k:k), 3,                            &
!            (/ 'VelocityX', 'VelocityY', 'VelocityZ' /) )
!#endif
        end if
#endif
    end do
    deallocate(ui,vi,wi)
else
    write(*,*) 'Error: itype not specified properly to inst_write!'
    stop
end if

deallocate(w_uv)
nullify(x,y,z,zw)

#ifdef PPLVLSET
contains
!*******************************************************************************
subroutine force_tot()
!*******************************************************************************
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
#endif
implicit none

! Zero bogus values
fx(:,:,nz) = 0._rprec
fy(:,:,nz) = 0._rprec
fz(:,:,nz) = 0._rprec

!  Sum both the induced and applied forces
allocate(fx_tot(nx,ny,nz), fy_tot(nx,ny,nz), fz_tot(nx,ny,nz))

#ifdef PPTURBINES
fx_tot = fxa(1:nx,1:ny,1:nz)
fy_tot = fya(1:nx,1:ny,1:nz)
fz_tot = fza(1:nx,1:ny,1:nz)

#elif defined(PPATM)
fx_tot = fxa(1:nx,1:ny,1:nz)
fy_tot = fya(1:nx,1:ny,1:nz)
fz_tot = fza(1:nx,1:ny,1:nz)

#elif defined(PPLVLSET)
fx_tot = fx(1:nx,1:ny,1:nz)+fxa(1:nx,1:ny,1:nz)
fy_tot = fy(1:nx,1:ny,1:nz)+fya(1:nx,1:ny,1:nz)
fz_tot = fz(1:nx,1:ny,1:nz)+fza(1:nx,1:ny,1:nz)
#else
fx_tot = 0._rprec
fy_tot = 0._rprec
fz_tot = 0._rprec
#endif

#ifdef PPMPI
!  Sync forces
call mpi_sync_real_array( fx_tot, 1, MPI_SYNC_DOWN )
call mpi_sync_real_array( fy_tot, 1, MPI_SYNC_DOWN )
call mpi_sync_real_array( fz_tot, 1, MPI_SYNC_DOWN )
#endif

! Put fz_tot on uv-grid
fz_tot(1:nx,1:ny,1:nz) = interp_to_uv_grid( fz_tot(1:nx,1:ny,1:nz), 1 )

return
end subroutine force_tot
#endif

!*******************************************************************************
!subroutine pressure_sync()
!!*******************************************************************************
!use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
!use param, only : ld
!implicit none
!
!! Reset bogus values
!p(:,:,nz) = p(:,:,nz-1)
!dpdx(:,:,nz) = dpdx(:,:,nz-1)
!dpdy(:,:,nz) = dpdy(:,:,nz-1)
!dpdz(:,:,nz) = dpdz(:,:,nz-1)
!
!#ifdef PPMPI
!!  Sync pressure
!call mpi_sync_real_array( p, 0 , MPI_SYNC_DOWN )
!call mpi_sync_real_array( dpdx, 1 , MPI_SYNC_DOWN )
!call mpi_sync_real_array( dpdy, 1 , MPI_SYNC_DOWN )
!call mpi_sync_real_array( dpdz, 1 , MPI_SYNC_DOWN )
!#endif
!
!return
!end subroutine pressure_sync
!
!!*******************************************************************************
!subroutine RHS_sync()
!!*******************************************************************************
!use param, only : ld
!use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
!implicit none
!
!! Reset bogus values
!RHSx(:,:,nz) = RHSx(:,:,nz-1)
!RHSy(:,:,nz) = RHSy(:,:,nz-1)
!RHSz(:,:,nz) = RHSz(:,:,nz-1)
!
!#ifdef PPMPI
!!  Sync RHS
!call mpi_sync_real_array( RHSx, 0 , MPI_SYNC_DOWN )
!call mpi_sync_real_array( RHSy, 0 , MPI_SYNC_DOWN )
!call mpi_sync_real_array( RHSz, 0 , MPI_SYNC_DOWN )
!#endif
!
!return
!end subroutine RHS_sync

end subroutine inst_write

!*******************************************************************************
subroutine checkpoint ()
!*******************************************************************************
use iwmles, only : iwm_checkPoint
use param, only : nz, checkpoint_file, tavg_calc, lbc_mom, L_x, L_y, L_z, path
#ifdef PPMPI
use param, only : comm, ierr
#endif
use sim_param, only : u, v, w, RHSx, RHSy, RHSz
use sgs_param, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN, lagran_dt
use param, only : jt_total, total_time, total_time_dim, dt,                    &
    use_cfl_dt, cfl, write_endian, inflow_type
use cfl_util, only : get_max_cfl
use string_util, only : string_concat
#ifdef PPTURBINES
use turbines, only : turbines_checkpoint
#endif
#ifdef PPSCALARS
use scalars, only : scalars_checkpoint
#endif
#ifdef PPATM
use atm_lesgo_interface, only : atm_lesgo_checkpoint
#endif
use coriolis, only : coriolis_finalize

! HIT Inflow
#ifdef PPHIT
use hit_inflow, only : hit_write_restart
#endif

implicit none
character(64) :: fname
real(rprec) :: cfl_w

#if defined(PPLES_GPU)
! GPU build: the restart state (u/v/w, AB2 RHS history, Lagrangian SGS
! accumulators) lives on the device; the host copies were last written at
! initialization. Refresh them so vel.out carries the evolved state.
! Infrequent path (every checkpoint_nskip steps + end of run).
!$acc wait(1)
!$acc update self(u, v, w, RHSx, RHSy, RHSz, Cs_opt2, F_LM, F_MM, F_QN, F_NN)
#endif

fname = checkpoint_file
#ifdef PPMPI
call string_concat( fname, '.c', coord )
#endif

!  Open vel.out (lun_default in io) for final output
open(11, file=fname, form='unformatted', convert=write_endian,                 &
    status='unknown', position='rewind')
write (11) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),                        &
    RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),                      &
    Cs_opt2(:,:,1:nz), F_LM(:,:,1:nz), F_MM(:,:,1:nz),                         &
    F_QN(:,:,1:nz), F_NN(:,:,1:nz)
close(11)

! Open grid.out for final output
if (coord == 0) then
    open(11, file= path // 'grid.out', form='unformatted', convert=write_endian)
    write(11) nproc, Nx, Ny, Nz, L_x, L_y, L_z
    close(11)
end if

#ifdef PPMPI
call mpi_barrier( comm, ierr )
#endif

! Checkpoint time averaging restart data
if ( tavg_calc .and. tavg%initialized ) call tavg%checkpoint()

! Write time and current simulation state
! cfl_w is the serialized CFL value: either the configured CFL target or the
! current maximum CFL when the run advances with a fixed dt.
if( use_cfl_dt ) then
    cfl_w = cfl
else
    cfl_w = get_max_cfl()
end if

! Checkpoint integral-wall-model state when that lower boundary model is active.
if(lbc_mom==3)then
    if (coord == 0) call iwm_checkPoint()
end if

#ifdef PPHIT
    if (coord == 0 .and. inflow_type == 2) call hit_write_restart()
#endif

#ifdef PPTURBINES
call turbines_checkpoint
#endif

#ifdef PPSCALARS
call scalars_checkpoint
#endif

#ifdef PPATM
! Keep ATM controller, induced-velocity, geometry, and structural histories at
! the same logical timestep as vel.out and total_time.dat.
call atm_lesgo_checkpoint(jt_total, total_time)
#endif

call coriolis_finalize()

!  Update total_time.dat after simulation
if (coord == 0) then
    !--only do this for true final output, not intermediate recording
    open (1, file=fcumulative_time)
    write(1, '(i0,1x,5(es26.17e3,1x))') jt_total, total_time, total_time_dim,   &
        dt, cfl_w, lagran_dt
    close(1)
end if

end subroutine checkpoint

!*******************************************************************************
subroutine output_final()
!*******************************************************************************
use param, only : tavg_calc
implicit none

! Perform final checkpoing
call checkpoint()

!  Check if average quantities are to be recorded
if (tavg_calc .and. tavg%initialized ) call tavg%finalize()

end subroutine output_final

!*******************************************************************************
subroutine output_init ()
!*******************************************************************************
!
!  This subroutine allocates the memory for arrays used for statistical
!  calculations
!
use param, only : dx, dy, dz, nz, lbz
use param, only : point_calc, point_nloc, point_loc
use param, only : xplane_calc, xplane_nloc, xplane_loc
use param, only : yplane_calc, yplane_nloc, yplane_loc
use param, only : zplane_calc, zplane_nloc, zplane_loc
use grid_m, only : grid
use functions, only : cell_indx
use stat_defs, only : point, xplane, yplane, zplane
implicit none

integer :: i,j,k
real(rprec), pointer, dimension(:) :: x,y,z


#ifdef PPMPI
! This adds one more element to the last processor (which contains an extra one)
! Processor nproc-1 has data from 1:nz
! Rest of processors have data from 1:nz-1
if ( coord == nproc-1 ) then
    nz_end = 0
else
    nz_end = 1
end if
#else
nz_end = 0
#endif

nullify(x,y,z)

x => grid % x
y => grid % y
z => grid % z

! Initialize information for x-planar stats/data
if (xplane_calc) then
    allocate(xplane(xplane_nloc))
    xplane(:) % istart = -1
    xplane(:) % ldiff = 0.

    !  Compute istart and ldiff
    do i = 1, xplane_nloc
        xplane(i) % istart = cell_indx('i', dx, xplane_loc(i))
        xplane(i) % ldiff = xplane_loc(i) - x(xplane(i) % istart)
    end do
end if

! Initialize information for y-planar stats/data
if (yplane_calc) then
    allocate(yplane(yplane_nloc))
    yplane(:) % istart = -1
    yplane(:) % ldiff = 0.

    !  Compute istart and ldiff
    do j = 1, yplane_nloc
        yplane(j) % istart = cell_indx('j', dy, yplane_loc(j))
        yplane(j) % ldiff = yplane_loc(j) - y(yplane(j) % istart)
    end do
end if

! Initialize information for z-planar stats/data
if(zplane_calc) then
    allocate(zplane(zplane_nloc))

    !  Initialize
    zplane(:) % istart = -1
    zplane(:) % ldiff = 0.
    zplane(:) % coord = -1

    !  Compute istart and ldiff
    do k = 1, zplane_nloc

#ifdef PPMPI
        if (zplane_loc(k) >= z(1) .and. zplane_loc(k) < z(nz)) then
            zplane(k) % coord = coord
            zplane(k) % istart = cell_indx('k',dz,zplane_loc(k))
            zplane(k) % ldiff = zplane_loc(k) - z(zplane(k) % istart)
        end if
#else
        zplane(k) % coord = 0
        zplane(k) % istart = cell_indx('k',dz,zplane_loc(k))
        zplane(k) % ldiff = zplane_loc(k) - z(zplane(k) % istart)
#endif
    end do
end if

!  Open files for instantaneous writing
if (point_calc) then
    allocate(point(point_nloc))

    !  Intialize the coord values
    ! (-1 shouldn't be used as coord so initialize to this)
    point % coord=-1
    point % fid = -1

    do i = 1, point_nloc
        !  Find the processor in which this point lives
#ifdef PPMPI
        if (point_loc(i)%xyz(3) >= z(1) .and. point_loc(i)%xyz(3) < z(nz)) then
#endif

            point(i) % coord = coord

            point(i) % istart = cell_indx('i',dx,point_loc(i)%xyz(1))
            point(i) % jstart = cell_indx('j',dy,point_loc(i)%xyz(2))
            point(i) % kstart = cell_indx('k',dz,point_loc(i)%xyz(3))

            point(i) % xdiff = point_loc(i)%xyz(1) - x(point(i) % istart)
            point(i) % ydiff = point_loc(i)%xyz(2) - y(point(i) % jstart)
            point(i) % zdiff = point_loc(i)%xyz(3) - z(point(i) % kstart)

#ifdef PPMPI
        end if
#endif
    end do
end if

nullify(x,y,z)

end subroutine output_init

end module io
