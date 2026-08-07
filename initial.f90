!!
!!  Copyright (C) 2009-2016  Johns Hopkins University
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
!!

!*******************************************************************************
subroutine initial()
!*******************************************************************************
use iwmles, only : iwm_read_checkPoint
use types,only:rprec
use param, only : BOGUS, L_x, L_y, L_z, checkpoint_file, coord,               &
                  cumulative_time, inflow_type, inflow_velocity, inilag,      &
                  initu, lbc_mom, lbz, nproc, nx, ny, nz, path, read_endian
use sim_param, only : u, v, w, RHSx, RHSy, RHSz
use sgs_param, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN
#ifdef PPDYN_TN
use sgs_param, only : F_ee2, F_deedt2, ee_past
#endif
#if (defined(PPTURBINES) && !defined(PPATM)) || defined(PPLVLSET)
use sim_param, only : fxa, fya, fza
#endif
#ifdef PPLVLSET
use sim_param, only : fx, fy, fz
use sgs_param, only : Beta,Tn_all
use level_set_base, only : phi
use messages, only : error
#endif
use string_util, only : string_concat
#ifdef PPMPI
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWNUP
#endif
#ifdef PPSCALARS
use scalars, only : ic_scal, theta
#endif

implicit none

character (64) :: fname

#ifdef PPDYN_TN
logical :: exst
character (64) :: fname_dyn_tn
#endif

integer::jz

! Flag to identify if file exists
logical :: file_flag
logical :: interp_flag
logical :: iwm_file_flag ! integral-wall-model restart file flag
#ifdef PPLVLSET
integer, parameter :: lvlset_restart_magic=20260807
integer, parameter :: lvlset_restart_version=1
character(64) :: fname_lvlset_restart
integer :: restart_magic,restart_version,restart_nx,restart_ny,restart_nz
integer :: restart_nproc
real(rprec) :: restart_phi_signature(4),current_phi_signature(4)
logical :: lvlset_restart_exists
#endif

#if defined(PPTURBINES) && !defined(PPATM)
fxa = 0._rprec
fya = 0._rprec
fza = 0._rprec
#endif

#ifdef PPLVLSET
fx = 0._rprec; fy = 0._rprec; fz = 0._rprec
fxa = 0._rprec;  fya = 0._rprec; fza = 0._rprec
#endif

#ifdef PPDYN_TN
! Will be over-written if read from dyn_tn.out files
ee_past = 0.1_rprec; F_ee2 = 10.0_rprec; F_deedt2 = 10000.0_rprec
fname_dyn_tn = path // 'dyn_tn.out'
#ifdef PPMPI
call string_concat( fname_dyn_tn, '.c', coord )
#endif
#endif

fname = checkpoint_file

#ifdef PPMPI
call string_concat( fname, '.c', coord )
#endif

! check iwm restart file
inquire (file='iwm_checkPoint.dat', exist=iwm_file_flag)
if (iwm_file_flag) then
    if (lbc_mom == 3) then
        if (coord==0) call iwm_read_checkPoint()
    end if
end if

! Check if file exists and read from it
! if not then create new IC
inquire( file=fname, exist=file_flag )
interp_flag = check_for_interp()
if (file_flag .and. .not.interp_flag) then
    initu = .true.
    inilag = .false.
else
    initu = .false.
    inilag = .true.
end if

if (initu) then
    if (coord == 0) write(*,*) '--> Reading initial velocity field from file'
    call ic_file
elseif (interp_flag) then
    if (coord == 0) write(*,*) '--> Interpolating initial velocity field from file'
    call ic_interp
else if (inflow_type == 1) then
    if (coord == 0) write(*,*) '--> Creating initial uniform velocity field'
    call ic_uniform
#ifdef PPSCALARS
    ! Scalar initial conditions are handled by ic_scal below.  A passive or
    ! prescribed scalar profile can be advanced with uniform velocity inflow.
    if (coord == 0) write(*,*) '--> Initializing scalar field with uniform velocity inflow'
#endif
else if (lbc_mom==1) then
    if (coord == 0) write(*,*) '--> Creating initial laminar profile ',&
        'field with DNS BCs'
    call ic_dns()
else
    if (coord == 0) write(*,*) '--> Creating initial boundary layer velocity ',&
    'field with LES BCs'
    call ic_les()
end if

#ifdef PPLVLSET
if (initu .and. allocated(Beta)) then
    fname_lvlset_restart=path//'lvlset_sgs_restart.out'
#ifdef PPMPI
    call string_concat(fname_lvlset_restart,'.c',coord)
#endif
    inquire(file=fname_lvlset_restart,exist=lvlset_restart_exists)
    if (.not. lvlset_restart_exists) call error('initial',                   &
        'missing Level Set model-4/5 restart state: '//trim(fname_lvlset_restart))
    open(14,file=fname_lvlset_restart,form='unformatted',convert=read_endian,&
        status='old',action='read')
    read(14) restart_magic,restart_version,restart_nx,restart_ny,restart_nz,&
        restart_nproc,restart_phi_signature
    if (restart_magic /= lvlset_restart_magic .or.                         &
        restart_version /= lvlset_restart_version) call error('initial',   &
        'unsupported Level Set SGS restart format')
    if (restart_nx /= nx .or. restart_ny /= ny .or. restart_nz /= nz .or. &
        restart_nproc /= nproc) call error('initial',                      &
        'Level Set SGS restart grid/decomposition mismatch')
    current_phi_signature=(/sum(phi(1:nx,:,1:nz-1)),                       &
        sum(abs(phi(1:nx,:,1:nz-1))),minval(phi(1:nx,:,1:nz-1)),          &
        maxval(phi(1:nx,:,1:nz-1))/)
    if (any(abs(current_phi_signature-restart_phi_signature) >             &
        1000._rprec*epsilon(1._rprec)*                                    &
        max(1._rprec,abs(restart_phi_signature)))) call error('initial',   &
        'Level Set SGS restart geometry signature mismatch')
    read(14) Beta(:,:,1:nz),Tn_all(:,:,1:nz)
    close(14)
end if
#endif

#ifdef PPDYN_TN
! Read dynamic timescale running averages from file
if (cumulative_time) then
    inquire (file=fname_dyn_tn, exist=exst)
    if (exst) then
        open(13,file=fname_dyn_tn,form='unformatted', convert=read_endian)
        read(13) F_ee2(:,:,1:nz), F_deedt2(:,:,1:nz), ee_past(:,:,1:nz)
    else
        write(*,*) trim(fname_dyn_tn), ' not found - using default values'
    end if
    close(13)
end if
#endif

#ifdef PPSCALARS
call ic_scal(interp_flag)
#endif

! call mpi_barrier(comm, ierr)
! stop

! Write averaged vertical profiles to standard output
do jz = 1, nz
#ifdef PPSCALARS
    write(6,7780) jz+coord*(nz-1), sum (u(1:nx, :, jz)) / (nx * ny),               &
                  sum (v(1:nx, :, jz)) / (nx * ny),                            &
                  sum (w(1:nx, :, jz)) / (nx * ny),                            &
                  sum (theta(1:nx, :, jz)) / (nx * ny)
#else
    write(6,7780) jz+coord*(nz-1), sum (u(1:nx, :, jz)) / (nx * ny),               &
                  sum (v(1:nx, :, jz)) / (nx * ny),                            &
                  sum (w(1:nx, :, jz)) / (nx * ny)
#endif
end do
#ifdef PPSCALARS
7780 format('jz, ubar, vbar, wbar, thetabar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4,1x,F9.4))
#else
7780 format('jz, ubar, vbar, wbar:',(1x,I3,1x,F9.4,1x,F9.4,1x,F9.4))
#endif

#ifdef PPMPI
! Exchange ghost node information for u, v, and w
call mpi_sync_real_array( u, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( v, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( w, 0, MPI_SYNC_DOWNUP )

! Restart files store physical planes 1:nz, not the lower overlap plane at
! index 0. Reconstruct every saved history used across a rank interface before
! the first resumed AB2 or Lagrangian SGS update.
if (initu) then
    call mpi_sync_real_array( RHSx, 0, MPI_SYNC_DOWNUP )
    call mpi_sync_real_array( RHSy, 0, MPI_SYNC_DOWNUP )
    call mpi_sync_real_array( RHSz, 0, MPI_SYNC_DOWNUP )
    call mpi_sync_real_array( F_LM, 0, MPI_SYNC_DOWNUP )
    call mpi_sync_real_array( F_MM, 0, MPI_SYNC_DOWNUP )
    call mpi_sync_real_array( F_QN, 0, MPI_SYNC_DOWNUP )
    call mpi_sync_real_array( F_NN, 0, MPI_SYNC_DOWNUP )
endif

!--set 0-level velocities to BOGUS
if (coord == 0) then
    u(:, :, lbz) = BOGUS
    v(:, :, lbz) = BOGUS
    w(:, :, lbz) = BOGUS
end if
#endif


contains


!*******************************************************************************
subroutine ic_uniform()
!*******************************************************************************
! This subroutine creates a uniform initial condition without turbulence.
!
implicit none

u = inflow_velocity
v = 0._rprec
w = 0._rprec

end subroutine ic_uniform

!*******************************************************************************
function check_for_interp() result(flag)
!*******************************************************************************
use param, only : path
integer :: Nx_f, Ny_f, Nz_f, nproc_f
real(rprec) :: Lx_f, Ly_f, Lz_f
logical :: exst, flag

inquire (file=path // 'grid.out', exist=exst)
if (exst) then
    open(12, file= path // 'grid.out', form='unformatted', convert=read_endian)
    read(12) nproc_f, Nx_f, Ny_f, Nz_f, Lx_f, Ly_f, Lz_f
    close(12)
    if (nproc_f == nproc .and. Nx_f == Nx .and. Ny_f == Ny .and. Nz_f == Nz    &
        .and. Lx_f == L_x .and. Ly_f == L_y .and. Lz_f == L_z) then
        flag = .false.
    else
        flag = .true.
    end if
else
    flag = .false.
end if

end function check_for_interp

!*******************************************************************************
subroutine ic_file()
!*******************************************************************************
! This subroutine reads the initial conditions from the checkpoint file.
!
open(12, file=fname, form='unformatted', convert=read_endian)

read(12) u(:, :, 1:nz), v(:, :, 1:nz), w(:, :, 1:nz),                          &
         RHSx(:, :, 1:nz), RHSy(:, :, 1:nz), RHSz(:, :, 1:nz),                 &
         Cs_opt2(:,:,1:nz), F_LM(:,:,1:nz), F_MM(:,:,1:nz),                    &
         F_QN(:,:,1:nz), F_NN(:,:,1:nz)

close(12)

end subroutine ic_file

!*******************************************************************************
subroutine ic_interp()
!*******************************************************************************
! This subroutine reads the initial conditions from a checkpoint file and
! interpolates onto the current grid
!
use param, only : path
use grid_m, only : grid
use functions, only : binary_search
integer :: nproc_f, Nx_f, Ny_f, Nz_f
real(rprec) :: Lx_f, Ly_f, Lz_f
integer :: i, j, k, z1, z2, ld_f, lh_f, Nz_tot_f
real(rprec) :: dx_f, dy_f, dz_f
integer :: i1, i2, j1, j2, k1, k2
real(rprec) :: ax, ay, az, bx, by, bz, xx, yy
real(rprec), allocatable, dimension(:) :: x_f, y_f!, z_f, zw_f
real(rprec), allocatable, dimension(:,:,:) :: u_f, v_f, w_f
character(64) :: ff
integer :: npr1, npr2, nproc_r, Nz_tot_r
real(rprec), allocatable, dimension(:) :: z_r, zw_r
! real(rprec), allocatable, dimension(:,:,:) :: u_r, v_r, w_r

! Read grid information from file
open(12, file=path // 'grid.out', form='unformatted', convert=read_endian)
read(12) nproc_f, Nx_f, Ny_f, Nz_f, Lx_f, Ly_f, Lz_f
close(12)

! Compute intermediate values
lh_f = nx_f/2+1
ld_f = 2*lh_f
Nz_tot_f = ( nz_f - 1 ) * nproc_f + 1
dx_f = Lx_f / nx_f
dy_f = Ly_f / ny_f
dz_f = Lz_f / (nz_tot_f - 1)

! Figure out which processors actually need to be read
npr1 = max(floor(grid%z(lbz)/dz_f/Nz_f), 0)
npr2 = min(ceiling(grid%z(nz)/dz_f/Nz_f), nproc_f-1)
nproc_r = npr2-npr1+1
Nz_tot_r = ( nz_f - 1 ) * nproc_r + 1
! write(*,*) coord, npr1, npr2, nproc_r, Nz_tot_r

! Create file grid
allocate( z_r(nz_tot_r), zw_r(nz_tot_r))
do i = 1, nz_tot_r
    zw_r(i) = (i - 1 + npr1*(nz_f-1)) * dz_f
    z_r(i) = zw_r(i) + 0.5*dz_f
end do

! Create file grid
allocate( x_f(nx_f), y_f(ny_f))
do i = 1, nx_f
    x_f(i) = (i-1) * Lx_f/(nx_f)
end do
do i = 1, ny_f
    y_f(i) = (i-1) * Ly_f/ny_f
end do

! Read velocities from file
allocate( u_f(ld_f, ny_f, nz_tot_r) )
allocate( v_f(ld_f, ny_f, nz_tot_r) )
allocate( w_f(ld_f, ny_f, nz_tot_r) )

! Loop through all levels
do i = 1, nproc_r
    ! Set level bounds
    z1 = nz_tot_f / nproc_f * (i-1) + 1
    z2 = nz_tot_f / nproc_f * i  + 1

    ! Read from file
    write(ff,*) i+npr1-1
    ff = path // "vel.out.c"//trim(adjustl(ff))
    open(12, file=ff,  action='read', form='unformatted')
    read(12) u_f(:, :, z1:z2), v_f(:, :, z1:z2), w_f(:, :, z1:z2)
    close(12)
end do

! Calculate velocities
do i = 1, nx
    xx = grid%x(i) - floor(grid%x(i)/Lx_f) * Lx_f
    i1 = binary_search(x_f, xx)
    i2 = mod(i1, nx_f) + 1
    bx = (xx - x_f(i1)) / dx_f
    ax = 1._rprec - bx
    do j = 1, ny
        yy = grid%y(j) - floor(grid%y(j)/Ly_f) * Ly_f
        j1 = binary_search(y_f, yy)
        j2 = mod(j1, ny_f) + 1
        by = (yy - y_f(j1)) / dy_f
        ay = 1._rprec - by
        do k = 1, nz
            ! for u and v
            k1 = binary_search(z_r, grid%z(k))
            k2 = k1 + 1
            if (k1 == nz_tot_r) then
                u(i,j,k) = ax*ay*u_f(i1,j1,k1) + bx*ay*u_f(i2,j1,k1)           &
                         + ax*by*u_f(i1,j2,k1) + bx*by*u_f(i2,j2,k1)
                v(i,j,k) = ax*ay*v_f(i1,j1,k1) + bx*ay*v_f(i2,j1,k1)           &
                         + ax*by*v_f(i1,j2,k1) + bx*by*v_f(i2,j2,k1)
            else if (k1 == 0) then
                u(i,j,k) = ax*ay*u_f(i1,j1,k2) + bx*ay*u_f(i2,j1,k2)           &
                         + ax*by*u_f(i1,j2,k2) + bx*by*u_f(i2,j2,k2)
                v(i,j,k) = ax*ay*v_f(i1,j1,k2) + bx*ay*v_f(i2,j1,k2)           &
                         + ax*by*v_f(i1,j2,k2) + bx*by*v_f(i2,j2,k2)
            else
                bz = (grid%z(k) - z_r(k1)) / dz_f
                az = 1._rprec - bz
                u(i,j,k) = ax*ay*az*u_f(i1,j1,k1) + bx*ay*az*u_f(i2,j1,k1)     &
                         + ax*by*az*u_f(i1,j2,k1) + bx*by*az*u_f(i2,j2,k1)     &
                         + ax*ay*bz*u_f(i1,j1,k2) + bx*ay*bz*u_f(i2,j1,k2)     &
                         + ax*by*bz*u_f(i1,j2,k2) + bx*by*bz*u_f(i2,j2,k2)
                v(i,j,k) = ax*ay*az*v_f(i1,j1,k1) + bx*ay*az*v_f(i2,j1,k1)     &
                         + ax*by*az*v_f(i1,j2,k1) + bx*by*az*v_f(i2,j2,k1)     &
                         + ax*ay*bz*v_f(i1,j1,k2) + bx*ay*bz*v_f(i2,j1,k2)     &
                         + ax*by*bz*v_f(i1,j2,k2) + bx*by*bz*v_f(i2,j2,k2)
            end if

            ! for w
            k1 = binary_search(zw_r, grid%zw(k))
            k2 = k1 + 1
            if (k1 == nz_tot_r) then
                w(i,j,k) = ax*ay*w_f(i1,j1,k1) + bx*ay*w_f(i2,j1,k1)           &
                         + ax*by*w_f(i1,j2,k1) + bx*by*w_f(i2,j2,k1)
            else if (k1 == 0) then
                w(i,j,k) = ax*ay*w_f(i1,j1,k2) + bx*ay*w_f(i2,j1,k2)           &
                         + ax*by*w_f(i1,j2,k2) + bx*by*w_f(i2,j2,k2)
            else
                bz = (grid%zw(k) - zw_r(k1)) / dz_f
                az = 1._rprec - bz
                w(i,j,k) = ax*ay*az*w_f(i1,j1,k1) + bx*ay*az*w_f(i2,j1,k1)     &
                         + ax*by*az*w_f(i1,j2,k1) + bx*by*az*w_f(i2,j2,k1)     &
                         + ax*ay*bz*w_f(i1,j1,k2) + bx*ay*bz*w_f(i2,j1,k2)     &
                         + ax*by*bz*w_f(i1,j2,k2) + bx*by*bz*w_f(i2,j2,k2)
            end if

        end do
    end do
end do

end subroutine ic_interp

!*******************************************************************************
subroutine ic_dns()
!*******************************************************************************
! This subroutine produces an initial condition for the boundary layer when
! using DNS boundary conditions.
!
use types,only:rprec
use param, only : coord, dz, L_z, nu_molec, nx, ny, nz, ubc_mom, ubot,        &
                  u_star, utop, z_i
use sim_param,only:u,v,w
implicit none

real(rprec), dimension(nz) :: ubar
real(rprec) :: rms, temp, z
integer :: jx, jy, jz
real(rprec) :: dummy_rand

! Calculate the average streamwise velocity based on height of first uvp point
! in wall units


if ( abs(ubot) > 0 .or. abs(utop) > 0 ) then
    !! linear laminar profile (couette) z and ubar are non-dimensional
    do jz = 1, nz
#ifdef PPMPI
        z = (coord*(nz-1) + real(jz,rprec) - 0.5_rprec) * dz
#else
        z = (real(jz,rprec) - 0.5_rprec) * dz ! non-dimensional
#endif
        ubar(jz)= (utop-ubot)/2*(2*z/L_z-1)**5 + (utop+ubot)/2
    end do
else
    !! parabolic laminar profile (channel) z and ubar are non-dimensional
    do jz = 1, nz
#ifdef PPMPI
        z = (coord*(nz-1) + real(jz,rprec) - 0.5_rprec) * dz
#else
        z = (real(jz,rprec) - 0.5_rprec) * dz
#endif
        ubar(jz)=(u_star*z_i/nu_molec) * z * (1._rprec - 0.5_rprec*z)
    end do
endif

! Get random seeds to populate the initial condition with noise
call init_random_seed

! Add noise to the velocity field
! the "default" rms of a unif variable is 0.289
rms = 0.2_rprec
do jz = 1, nz
    do jy = 1, ny
        do jx = 1, nx
            call random_number(dummy_rand)
            u(jx,jy,jz)=ubar(jz)+(rms/.289_rprec)*(dummy_rand-.5_rprec)/u_star
            call random_number(dummy_rand)
            v(jx,jy,jz) = 0._rprec+(rms/.289_rprec)*(dummy_rand-.5_rprec)/u_star
            call random_number(dummy_rand)
            w(jx,jy,jz) = 0._rprec+(rms/.289_rprec)*(dummy_rand-.5_rprec)/u_star
        end do
    end do
end do

! make sure w-mean is 0
temp=0._rprec
do jz = 1, nz
    do jy = 1, ny
        do jx = 1, nx
            temp = temp+w(jx,jy,jz)
        end do
    end do
end do
temp = temp/(nx*ny*nz)

do jz = 1, nz
   do jy = 1, ny
      do jx = 1, nx
         w(jx,jy,jz) = w(jx,jy,jz)-temp
      end do
   end do
end do

! Make sure field satisfies boundary conditions
w(:,:,1) = 0._rprec
w(:,:,nz) = 0._rprec
if (ubc_mom == 0) then
   u(:,:,nz) = u(:,:,nz-1)
   v(:,:,nz) = v(:,:,nz-1)
endif

end subroutine ic_dns

!*******************************************************************************
subroutine ic_les()
!*******************************************************************************
! This subroutine produces an initial condition for the boundary layer.
! A log profile is used that flattens at z=z_i. Noise is added to
! promote the generation of turbulence
!
use types,only:rprec
use param, only : coord, dz, lbc_mom, mean_p_force_x, mean_p_force_y, nproc,  &
                  nx, ny, nz, ubc_mom, u_star, vonk, z_i, zo
use sim_param, only : u, v, w
use messages, only : error
use coriolis, only : coriolis_forcing, G, alpha, fc
#if defined(PPTURBINES) && !defined(PPATM)
use turbines, only : turbine_vel_init
#endif
#ifdef PPSCALARS
use scalars, only : ic_no_vel_noise_z
#endif

implicit none
integer :: jz, jz_abs
real(rprec), dimension(nz) :: ubar, vbar
real(rprec) :: rms, sigma_rv, arg, z, mean_p_force_mag
real(rprec) :: gamma_e, u_par, u_perp
real(rprec), parameter :: Km = 5._rprec

character(*), parameter :: sub_name = 'ic'

#if defined(PPTURBINES) && !defined(PPATM)
real(rprec) :: zo_turbines = 0._rprec
#endif

do jz = 1, nz
#ifdef PPMPI
    z = (coord*(nz-1) + jz - 0.5_rprec) * dz
#else
    z = (jz - 0.5_rprec) * dz
#endif

    ! For channel flow, choose closest wall
    if(lbc_mom  > 0 .and. ubc_mom > 0) z = min(z, dz*nproc*(nz-1) - z)
    ! For upside-down half-channel, choose upper wall
    if(lbc_mom == 0 .and. ubc_mom > 0) z = dz*nproc*(nz-1) - z

    ! IC in equilibrium with rough surface (rough dominates in effective zo)
    arg = (1._rprec/vonk)*log(z/zo)

#ifdef PPLVLSET
    ! Kludge to adjust magnitude of velocity profile
    ! Not critical - may delete
    arg = 0.357*arg
#endif

#if defined(PPTURBINES) && !defined(PPATM)
    call turbine_vel_init (zo_turbines)
    arg = (1._rprec/vonk)*log(z/zo_turbines)
#endif

    mean_p_force_mag = sqrt(mean_p_force_x**2 + mean_p_force_y**2)
    if (coriolis_forcing > 0) then
        if (arg < G) then
            ubar(jz) = arg*cos(alpha)
            vbar(jz) = arg*sin(alpha)
        else
            ubar(jz) = G*cos(alpha)
            vbar(jz) = G*sin(alpha)
        endif
    else
        if (mean_p_force_mag > 0._rprec) then
            ubar(jz) = arg*mean_p_force_x/mean_p_force_mag
            vbar(jz) = arg*mean_p_force_y/mean_p_force_mag
        else
            ubar(jz) = arg
            vbar(jz) = 0._rprec
        endif
    end if
end do

rms = 3._rprec
sigma_rv = 0.289_rprec

! Fill u, v, and w with uniformly distributed random numbers between 0 and 1
call init_random_seed
call random_number(u)
call random_number(v)
call random_number(w)

! Center random number about 0 and rescale
u = rms / sigma_rv * (u - 0.5_rprec)
v = rms / sigma_rv * (v - 0.5_rprec)
w = rms / sigma_rv * (w - 0.5_rprec)


! Rescale noise depending on distance from wall and mean log profile
! z is in meters
do jz = 1, nz
#ifdef PPMPI
    jz_abs = coord * (nz-1) + jz
    z = (coord * (nz-1) + jz - 0.5_rprec) * dz * z_i
#else
    jz_abs = jz
    z = (jz-.5_rprec) * dz * z_i
#endif

    ! For channel flow, choose closest wall
    if(lbc_mom  > 0 .and. ubc_mom > 0) z = min(z, dz*nproc*(nz-1)*z_i - z)
    ! For upside-down half-channel, choose upper wall
    if(lbc_mom == 0 .and. ubc_mom > 0) z = dz*nproc*(nz-1)*z_i - z

#ifdef PPSCALARS
    if (z > ic_no_vel_noise_z*z_i) then
        u(:,:,jz) = ubar(jz)
        v(:,:,jz) = vbar(jz)
        w(:,:,jz) = 0._rprec
    elseif (z <= z_i) then
#else
    if (z <= z_i) then
#endif
        u(:,:,jz) = u(:,:,jz) * (1._rprec-z / z_i) + ubar(jz)
        v(:,:,jz) = v(:,:,jz) * (1._rprec-z / z_i) + vbar(jz)
        w(:,:,jz) = w(:,:,jz) * (1._rprec-z / z_i)
    else
        u(:,:,jz) = u(:,:,jz) * 0.01_rprec + ubar(jz)
        v(:,:,jz) = v(:,:,jz) * 0.01_rprec + vbar(jz)
        w(:,:,jz) = w(:,:,jz) * 0.01_rprec
    end if
end do

! Bottom boundary conditions
if (coord == 0) then
    w(:, :, 1) = 0._rprec
#ifdef PPMPI
    u(:, :, 0) = 0._rprec
    v(:, :, 0) = 0._rprec
    w(:, :, 0) = 0._rprec
#endif
end if

! Set upper boundary condition as zero for u, v, and w
#ifdef PPMPI
if (coord == nproc-1) then
#endif
    w(1:nx, 1:ny, nz) = 0._rprec
    u(1:nx, 1:ny, nz) = 0._rprec
    v(1:nx, 1:ny, nz) = 0._rprec
#ifdef PPMPI
end if
#endif

end subroutine ic_les

end subroutine initial
