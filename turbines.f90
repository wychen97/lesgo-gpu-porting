!!
!!  Copyright (C) 2010-2016  Johns Hopkins University
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
module turbines
!*******************************************************************************
! This module contains all of the subroutines associated with drag-disk turbines
!
! Navigation map:
!   - input and placement: turbines_init, turbines_nodes, place_turbines
!   - GPU metadata path: turbines_acc_metadata_init/finalize/sync_device_field
!   - timestep forcing: turbines_forcing_acc and turbines_forcing
!   - output/restart: turbines_checkpoint and turbines_finalize
!   - dynamic controls: turbine_vel_init and read_control_files
!
! This module owns the actuator-disk turbine path.  The actuator-line/ATM path
! is owned by actuator_turbine_model.f90 and atm_lesgo_interface.f90.

use types, only : rprec
use param
use grid_m, only : grid
use messages, only : error, warn
use string_util, only : string_splice
use turbine_indicator, only : turb_ind_func_t
use turbines_gpu, only : turbines_cuda_enabled, turbines_interp_w_to_uv_gpu
use functions, only : count_lines
use stat_defs, only : wind_farm
#ifdef PPMPI
use mpi_defs, only : MPI_SYNC_DOWN, MPI_SYNC_DOWNUP, mpi_sync_real_array
#endif

implicit none

save
private

public :: turbines_init, turbines_forcing, turbine_vel_init, turbines_finalize, &
    turbines_checkpoint, turbines_acc_available, turbines_forcing_acc

character (*), parameter :: mod_name = 'turbines'

! The following values are read from the input file
! number of turbines in the x-direction
integer, public :: num_x
! number of turbines in the y-direction
integer, public :: num_y
! baseline diameter in meters
real(rprec), public :: dia_all
! baseline height in meters
real(rprec), public :: height_all
! baseline thickness in meters
real(rprec), public :: thk_all
! orientation of turbines
integer, public :: orientation
! stagger percentage from baseline
real(rprec), public :: stag_perc
! angle from upstream (CCW from above, -x dir is zero)
real(rprec), public :: theta1_all
! angle above horizontal
real(rprec), public :: theta2_all
! thrust coefficient (default 1.33)
real(rprec), public :: Ct_prime
! Read parameters from input_turbines/param.dat
logical, public :: read_param
! Dynamically change theta1 from input_turbines/theta1.dat
logical, public :: dyn_theta1
! Dynamically change theta2 from input_turbines/theta2.dat
logical, public :: dyn_theta2
! Dynamically change Ct_prime from input_turbines/Ct_prime.dat
logical, public :: dyn_Ct_prime
! Use ADM with rotation
logical, public :: use_rotation = .false.
! Tip speed ratio for ADM with rotation
real(rprec), public :: tip_speed_ratio = 7
! disk-avg time scale in seconds (default 600)
real(rprec), public :: T_avg_dim
! filter size as multiple of grid spacing
real(rprec), public :: alpha1
real(rprec), public :: alpha2
! indicator function only includes values above this threshold
real(rprec), public :: filter_cutoff
! Correct ADM for filtered indicator function
logical, public :: adm_correction
! Number of timesteps between the output
integer, public :: tbase

! The following are derived from the values above
integer :: nloc             ! total number of turbines
real(rprec) :: sx           ! spacing in the x-direction, multiple of diameter
real(rprec) :: sy           ! spacing in the y-direction

! Arrays for interpolating dynamic controls
real(rprec), dimension(:,:), allocatable :: theta1_arr
real(rprec), dimension(:), allocatable :: theta1_time
real(rprec), dimension(:,:), allocatable :: theta2_arr
real(rprec), dimension(:), allocatable :: theta2_time
real(rprec), dimension(:,:), allocatable :: Ct_prime_arr
real(rprec), dimension(:), allocatable :: Ct_prime_time

! Input files
character(:), allocatable :: input_folder
character(:), allocatable :: param_dat, theta1_dat, theta2_dat, Ct_prime_dat

! Output files
character(:), allocatable :: output_folder
character(:), allocatable :: vel_top_dat , u_d_T_dat
integer, dimension(:), allocatable :: forcing_fid

! epsilon used for disk velocity time-averaging
real(rprec) :: eps

! Commonly used indices
integer :: i, j, k, i2, j2, k2, l, s
integer :: k_start, k_end

#ifdef PPLES_GPU
logical, save :: turbines_acc_ready = .false.
integer, save :: acc_nloc = 0, acc_total_nodes = 0
integer, save :: acc_nloc_capacity = 0, acc_node_capacity = 0
integer, save :: acc_use_rotation = 0
real(rprec), save :: acc_inv_tip_speed_ratio = 0._rprec
integer, allocatable :: acc_offset(:), acc_count(:)
integer, allocatable :: acc_i(:), acc_j(:), acc_k(:)
integer, allocatable :: acc_icp(:), acc_jcp(:), acc_kcp(:), acc_center(:)
real(rprec), allocatable :: acc_ind(:), acc_ind_t(:)
real(rprec), allocatable :: acc_e1(:), acc_e2(:), acc_e3(:)
real(rprec), allocatable :: acc_nhat1(:), acc_nhat2(:), acc_nhat3(:)
real(rprec), allocatable :: acc_dia(:), acc_ct(:), acc_ud(:), acc_udt(:), acc_fn(:)
real(rprec), allocatable :: acc_disk(:), acc_uc(:), acc_vc(:), acc_wc(:)
real(rprec), allocatable :: acc_w_uv(:,:,:), acc_fza_uv(:,:,:)
!$acc declare create(acc_offset, acc_count, acc_i, acc_j, acc_k)
!$acc declare create(acc_icp, acc_jcp, acc_kcp, acc_center)
!$acc declare create(acc_ind, acc_ind_t, acc_e1, acc_e2, acc_e3)
!$acc declare create(acc_nhat1, acc_nhat2, acc_nhat3)
!$acc declare create(acc_dia, acc_ct, acc_ud, acc_udt, acc_fn)
!$acc declare create(acc_disk, acc_uc, acc_vc, acc_wc, acc_w_uv, acc_fza_uv)
!$acc declare create(acc_use_rotation, acc_inv_tip_speed_ratio)
!$acc declare create(acc_nloc, acc_total_nodes)
#endif

contains

!*******************************************************************************
subroutine turbines_init()
!*******************************************************************************
!
! This subroutine creates the 'turbine' folder and starts the turbine forcing
! output files. It also creates the indicator function (Gaussian-filtered from
! binary locations - in or out) and sets values for turbine type
! (node locations, etc)
!
implicit none

real(rprec), pointer, dimension(:) :: x,y,z
character (*), parameter :: sub_name = mod_name // '.turbines_init'
integer :: fid
real(rprec) :: T_avg_dim_file, delta1, delta2
logical :: test_logical, exst
character (100) :: string1

! Set pointers
nullify(x,y,z)
x => grid % x
y => grid % y
z => grid % z

! Input/Output file names
allocate(input_folder, source = 'input_turbines/')
allocate(param_dat, source = path // input_folder // 'param.dat')
allocate(theta1_dat, source = path // input_folder // 'theta1.dat')
allocate(theta2_dat, source = path // input_folder // 'theta2.dat')
allocate(Ct_prime_dat, source = path // input_folder // 'Ct_prime.dat')
allocate(output_folder, source = 'turbine/')
allocate(vel_top_dat, source = path // output_folder // 'vel_top.dat')
allocate(u_d_T_dat, source = path // output_folder // 'u_d_T.dat')

! Allocate and initialize
nloc = num_x*num_y
nullify(wind_farm%turbine)
allocate(wind_farm%turbine(nloc))
allocate(forcing_fid(nloc))

! Create turbine directory
call system('mkdir -vp ' // path // output_folder)

! Non-dimensionalize length values by z_i
height_all = height_all / z_i
dia_all = dia_all / z_i
thk_all = thk_all / z_i

! Spacing between turbines (as multiple of mean diameter)
sx = L_x / (num_x * dia_all )
sy = L_y / (num_y * dia_all )

! Place the turbines and specify some parameters
call place_turbines

! Resize thickness to capture at least on plane of gridpoints
! and set baseline values for size
do k = 1, nloc
    wind_farm%turbine(k)%thk = max(wind_farm%turbine(k)%thk, dx * 1.01)
    ! wind_farm%turbine(k)%vol_c = dx*dy*dz/(pi/4.*(wind_farm%turbine(k)%dia)**2 &
    !     * wind_farm%turbine(k)%thk)
end do

! Specify starting and ending indices for the processor
#ifdef PPMPI
k_start = 1+coord*(nz-1)
k_end = nz-1+coord*(nz-1)
#else
k_start = 1
k_end = nz
#endif

! Find the center of each turbine
do k = 1,nloc
    wind_farm%turbine(k)%icp = nint(wind_farm%turbine(k)%xloc/dx)
    wind_farm%turbine(k)%jcp = nint(wind_farm%turbine(k)%yloc/dy)
    wind_farm%turbine(k)%kcp = nint(wind_farm%turbine(k)%height/dz + 0.5)

    ! Check if turbine is the current processor
    test_logical = wind_farm%turbine(k)%kcp >= k_start .and.                   &
           wind_farm%turbine(k)%kcp<=k_end
    if (test_logical) then
        wind_farm%turbine(k)%center_in_proc = .true.
    else
        wind_farm%turbine(k)%center_in_proc = .false.
    end if

    ! Make kcp the local index
    wind_farm%turbine(k)%kcp = wind_farm%turbine(k)%kcp - k_start + 1

end do

! Read dynamic control input files
call read_control_files

!Compute a lookup table object for the indicator function
delta1 = alpha1 * sqrt(dx**2 + dy**2 + dz**2)
delta2 = alpha2 * sqrt(dx**2 + dy**2 + dz**2)
do s = 1, nloc
    call  wind_farm%turbine(s)%turb_ind_func%init(delta1, delta2,              &
            wind_farm%turbine(s)%thk, wind_farm%turbine(s)%dia)
end do

! Find turbine nodes - including filtered ind, n_hat, num_nodes, and nodes for
! each turbine. Each processor finds turbines in its domain
call turbines_nodes

! Read the time-averaged disk velocities from file if available
if (coord == 0) then
    inquire (file=u_d_T_dat, exist=exst)
    if (exst) then
        write(*,*) 'Reading from file ', trim(u_d_T_dat)
        open(newunit=fid, file=u_d_T_dat, status='unknown', form='formatted',  &
            position='rewind')
        do i=1,nloc
            read(fid,*) wind_farm%turbine(i)%u_d_T
        end do
        read(fid,*) T_avg_dim_file
        if (T_avg_dim_file /= T_avg_dim) then
            write(*,*) 'Time-averaging window does not match value in ',       &
                       trim(u_d_T_dat)
        end if
        close (fid)
    else
        write (*, *) 'File ', trim(u_d_T_dat), ' not found'
        write (*, *) 'Assuming u_d_T = -1. for all turbines'
        do k=1,nloc
            wind_farm%turbine(k)%u_d_T = -1.
        end do
    end if
end if

! Generate top of domain file
if (coord .eq. nproc-1) then
    open(newunit=fid, file=vel_top_dat, status='unknown', form='formatted',    &
        position='rewind')
    close(fid)
end if

! Generate the files for the turbine forcing output
if(coord==0) then
    do s=1,nloc
        call string_splice( string1, path // 'turbine/turbine_', s, '.dat' )
        open(newunit=forcing_fid(s), file=string1, status='unknown',           &
            form='formatted', position='append')
    end do
end if

#ifdef PPLES_GPU
call turbines_acc_metadata_init()
#endif

nullify(x,y,z)

end subroutine turbines_init

!*******************************************************************************
subroutine turbines_nodes
!*******************************************************************************
!
! This subroutine locates nodes for each turbine and builds the arrays: ind,
! n_hat, num_nodes, and nodes
!
use functions, only : cross_product
implicit none

character (*), parameter :: sub_name = mod_name // '.turbines_nodes'

real(rprec) :: rx,ry,rz,r,r_norm,r_disk

real(rprec), pointer :: p_xloc => null(), p_yloc => null(), p_height => null()
real(rprec), pointer :: p_dia => null(), p_thk => null()
real(rprec), pointer :: p_theta1 => null(), p_theta2 => null()
real(rprec), pointer :: p_nhat1 => null(), p_nhat2=> null(), p_nhat3 => null()
integer :: icp, jcp, kcp
integer :: imax, jmax, kmax
integer :: min_i, max_i, min_j, max_j, min_k, max_k
integer :: count_i, count_n
real(rprec), dimension(nz_tot) :: z_tot
real(rprec), dimension(3) :: temp_vec

#ifdef PPMPI
real(rprec), dimension(nloc) :: buffer_array
#endif
real(rprec), pointer, dimension(:) :: x, y, z

real(rprec) :: filt, filt_t, search_rad, filt_max
real(rprec), dimension(nloc) :: sumA, turbine_vol

nullify(x,y,z)

x => grid % x
y => grid % y
z => grid % z

sumA = 0._rprec

! z_tot for total domain (since z is local to the processor)
do k = 1,nz_tot
    z_tot(k) = (k - 0.5_rprec) * dz
end do

do s=1,nloc

    count_n = 0    !used for counting nodes for each turbine
    count_i = 1    !index count - used for writing to array "nodes"

    !set pointers
    p_xloc => wind_farm%turbine(s)%xloc
    p_yloc => wind_farm%turbine(s)%yloc
    p_height => wind_farm%turbine(s)%height
    p_dia => wind_farm%turbine(s)%dia
    p_thk => wind_farm%turbine(s)%thk
    p_theta1 => wind_farm%turbine(s)%theta1
    p_theta2 => wind_farm%turbine(s)%theta2
    p_nhat1 => wind_farm%turbine(s)%nhat(1)
    p_nhat2 => wind_farm%turbine(s)%nhat(2)
    p_nhat3 => wind_farm%turbine(s)%nhat(3)

    !identify "search area"
    search_rad = 0.5_rprec*p_dia + 3*max(alpha1, alpha2) * sqrt(dx**2 + dy**2 + dz**2)
    imax = min(int(search_rad/dx + 2), Nx/2)
    jmax = min(int(search_rad/dy + 2), Ny/2)
    kmax = int(search_rad/dz + 2)

    !determine unit normal vector for each turbine
    p_nhat1 = -cos(pi*p_theta1/180.)*cos(pi*p_theta2/180.)
    p_nhat2 = -sin(pi*p_theta1/180.)*cos(pi*p_theta2/180.)
    p_nhat3 = sin(pi*p_theta2/180.)

    !determine nearest (i,j,k) to turbine center
    icp = nint(p_xloc/dx)
    jcp = nint(p_yloc/dy)
    kcp = nint(p_height/dz + 0.5)

    !determine limits for checking i,j,k
    !due to spectral BCs, i and j may be < 1 or > nx,ny
    !the mod function accounts for this when these values are used
    min_i = icp-imax+1
    max_i = icp+imax
    min_j = jcp-jmax+1
    max_j = jcp+jmax
    min_k = max((kcp-kmax),1)
    max_k = min((kcp+kmax),nz_tot)
    wind_farm%turbine(s)%nodes_max(1) = min_i
    wind_farm%turbine(s)%nodes_max(2) = max_i
    wind_farm%turbine(s)%nodes_max(3) = min_j
    wind_farm%turbine(s)%nodes_max(4) = max_j
    wind_farm%turbine(s)%nodes_max(5) = min_k
    wind_farm%turbine(s)%nodes_max(6) = max_k

    ! check neighboring grid points
    ! update num_nodes, nodes, and ind for this turbine
    ! split domain between processors
    ! z(nz) and z(1) of neighboring coords match so each coord gets
    ! (local) 1 to nz-1
    wind_farm%turbine(s)%ind = 0._rprec
    wind_farm%turbine(s)%nodes = 0
    wind_farm%turbine(s)%num_nodes = 0
    count_n = 0
    count_i = 1

    ! Maximum value the filter takes (should be 1/volume)
    call wind_farm%turbine(s)%turb_ind_func%val(0._rprec, 0._rprec, filt_max, r_disk)

    do k=k_start,k_end  !global k
        do j=min_j,max_j
            do i=min_i,max_i
                ! vector from center point to this node is (rx,ry,rz)
                ! with length r
                if (i<1) then
                    i2 = mod(i+nx-1,nx)+1
                    rx = (x(i2)-L_x) - p_xloc
                elseif (i>nx) then
                    i2 = mod(i+nx-1,nx)+1
                    rx = (L_x+x(i2)) - p_xloc
                else
                    i2 = i
                    rx = x(i) - p_xloc
                end if
                if (j<1) then
                    j2 = mod(j+ny-1,ny)+1
                    ry = (y(j2)-L_y) - p_yloc
                elseif (j>ny) then
                    j2 = mod(j+ny-1,ny)+1
                    ry = (L_y+y(j2)) - p_yloc
                else
                    j2 = j
                    ry = y(j) - p_yloc
                end if
                rz = z_tot(k) - p_height
                r = sqrt(rx*rx + ry*ry + rz*rz)
                !length projected onto unit normal for this turbine
                r_norm = abs(rx*p_nhat1 + ry*p_nhat2 + rz*p_nhat3)
                !(remaining) length projected onto turbine disk
                r_disk = sqrt(r*r - r_norm*r_norm)
                ! get the filter value
                call wind_farm%turbine(s)%turb_ind_func%val(r_disk, r_norm, filt, filt_t)

                if ( filt > filter_cutoff * filt_max ) then
                    wind_farm%turbine(s)%ind(count_i) = filt
                    wind_farm%turbine(s)%ind_t(count_i) = filt_t
                    temp_vec(1) = rx-r_norm*p_nhat1
                    temp_vec(2) = ry-r_norm*p_nhat2
                    temp_vec(3) = rz-r_norm*p_nhat3
                    wind_farm%turbine(s)%e_theta(count_i,:) =                  &
                        cross_product(wind_farm%turbine(s)%nhat, temp_vec)
                    wind_farm%turbine(s)%e_theta(count_i,:) =                  &
                        wind_farm%turbine(s)%e_theta(count_i,:)                &
                        / sqrt(wind_farm%turbine(s)%e_theta(count_i,1)**2      &
                        + wind_farm%turbine(s)%e_theta(count_i,2)**2           &
                        + wind_farm%turbine(s)%e_theta(count_i,3)**2)
                    wind_farm%turbine(s)%nodes(count_i,1) = i2
                    wind_farm%turbine(s)%nodes(count_i,2) = j2
                    wind_farm%turbine(s)%nodes(count_i,3) = k-coord*(nz-1)!local
                    count_n = count_n + 1
                    count_i = count_i + 1
                    sumA(s) = sumA(s) + filt * dx * dy * dz
                end if
           end do
       end do
    end do
    wind_farm%turbine(s)%num_nodes = count_n

    ! Calculate turbine volume
    turbine_vol(s) = 0.25 * pi* p_dia**2 * p_thk

end do

! Sum the indicator function across all processors if using MPI
#ifdef PPMPI
buffer_array = sumA
call MPI_Allreduce(buffer_array, sumA, nloc, MPI_rprec, MPI_SUM, comm, ierr)
#endif

! Normalize the indicator function integrate to 1
do s = 1, nloc
    wind_farm%turbine(s)%ind=wind_farm%turbine(s)%ind(:)/sumA(s)
    wind_farm%turbine(s)%ind_t=wind_farm%turbine(s)%ind_t(:)/sumA(s)
end do

! Cleanup
nullify(x,y,z)

end subroutine turbines_nodes

#ifdef PPLES_GPU
!*******************************************************************************
logical function turbines_acc_available()
!*******************************************************************************
implicit none

turbines_acc_available = turbines_acc_ready

end function turbines_acc_available

!*******************************************************************************
subroutine turbines_acc_metadata_init()
!*******************************************************************************
implicit none

integer :: pos, node, needed_nloc, needed_total_nodes, needed_node_capacity
logical :: need_alloc

turbines_acc_ready = .false.

needed_nloc = nloc
needed_total_nodes = 0
do s = 1, nloc
    needed_total_nodes = needed_total_nodes + wind_farm%turbine(s)%num_nodes
end do
needed_node_capacity = max(1, needed_total_nodes)

need_alloc = .not. allocated(acc_offset) .or.                              &
    needed_nloc > acc_nloc_capacity .or.                                    &
    needed_node_capacity > acc_node_capacity .or.                            &
    .not. allocated(acc_w_uv) .or. .not. allocated(acc_fza_uv)

if (need_alloc) then
    call turbines_acc_finalize()
    acc_nloc_capacity = max(1, needed_nloc)
    acc_node_capacity = needed_node_capacity

    allocate(acc_offset(acc_nloc_capacity), acc_count(acc_nloc_capacity))
    allocate(acc_icp(acc_nloc_capacity), acc_jcp(acc_nloc_capacity),         &
        acc_kcp(acc_nloc_capacity), acc_center(acc_nloc_capacity))
    allocate(acc_nhat1(acc_nloc_capacity), acc_nhat2(acc_nloc_capacity),     &
        acc_nhat3(acc_nloc_capacity))
    allocate(acc_dia(acc_nloc_capacity), acc_ct(acc_nloc_capacity),          &
        acc_ud(acc_nloc_capacity), acc_udt(acc_nloc_capacity),               &
        acc_fn(acc_nloc_capacity))
    allocate(acc_disk(acc_nloc_capacity), acc_uc(acc_nloc_capacity),         &
        acc_vc(acc_nloc_capacity), acc_wc(acc_nloc_capacity))
    allocate(acc_i(acc_node_capacity), acc_j(acc_node_capacity),             &
        acc_k(acc_node_capacity))
    allocate(acc_ind(acc_node_capacity), acc_ind_t(acc_node_capacity))
    allocate(acc_e1(acc_node_capacity), acc_e2(acc_node_capacity),           &
        acc_e3(acc_node_capacity))
    allocate(acc_w_uv(ld, ny, lbz:nz), acc_fza_uv(ld, ny, lbz:nz))
end if

acc_nloc = needed_nloc
acc_total_nodes = needed_total_nodes

acc_offset = 0
acc_count = 0
acc_icp = 0
acc_jcp = 0
acc_kcp = 0
acc_center = 0
acc_nhat1 = 0._rprec
acc_nhat2 = 0._rprec
acc_nhat3 = 0._rprec
acc_dia = 0._rprec
acc_ct = 0._rprec
acc_ud = 0._rprec
acc_udt = 0._rprec
acc_fn = 0._rprec
acc_disk = 0._rprec
acc_uc = 0._rprec
acc_vc = 0._rprec
acc_wc = 0._rprec
acc_i = 0
acc_j = 0
acc_k = 0
acc_ind = 0._rprec
acc_ind_t = 0._rprec
acc_e1 = 0._rprec
acc_e2 = 0._rprec
acc_e3 = 0._rprec
if (need_alloc) then
    acc_w_uv = 0._rprec
    acc_fza_uv = 0._rprec
end if

pos = 1
do s = 1, acc_nloc
    acc_offset(s) = pos
    acc_count(s) = wind_farm%turbine(s)%num_nodes
    acc_icp(s) = wind_farm%turbine(s)%icp
    acc_jcp(s) = wind_farm%turbine(s)%jcp
    acc_kcp(s) = wind_farm%turbine(s)%kcp
    if (wind_farm%turbine(s)%center_in_proc) acc_center(s) = 1
    acc_nhat1(s) = wind_farm%turbine(s)%nhat(1)
    acc_nhat2(s) = wind_farm%turbine(s)%nhat(2)
    acc_nhat3(s) = wind_farm%turbine(s)%nhat(3)
    acc_dia(s) = wind_farm%turbine(s)%dia
    acc_ct(s) = wind_farm%turbine(s)%Ct_prime
    acc_ud(s) = wind_farm%turbine(s)%u_d
    acc_udt(s) = wind_farm%turbine(s)%u_d_T
    acc_fn(s) = wind_farm%turbine(s)%f_n
    do node = 1, acc_count(s)
        acc_i(pos) = wind_farm%turbine(s)%nodes(node,1)
        acc_j(pos) = wind_farm%turbine(s)%nodes(node,2)
        acc_k(pos) = wind_farm%turbine(s)%nodes(node,3)
        acc_ind(pos) = wind_farm%turbine(s)%ind(node)
        acc_ind_t(pos) = wind_farm%turbine(s)%ind_t(node)
        acc_e1(pos) = wind_farm%turbine(s)%e_theta(node,1)
        acc_e2(pos) = wind_farm%turbine(s)%e_theta(node,2)
        acc_e3(pos) = wind_farm%turbine(s)%e_theta(node,3)
        pos = pos + 1
    end do
end do

if (use_rotation) then
    acc_use_rotation = 1
    acc_inv_tip_speed_ratio = 1._rprec / tip_speed_ratio
else
    acc_use_rotation = 0
    acc_inv_tip_speed_ratio = 0._rprec
end if

if (acc_nloc > 0) then
    !$acc update device(acc_offset(1:acc_nloc), acc_count(1:acc_nloc))
    !$acc update device(acc_icp(1:acc_nloc), acc_jcp(1:acc_nloc),             &
    !$acc               acc_kcp(1:acc_nloc), acc_center(1:acc_nloc))
    !$acc update device(acc_nhat1(1:acc_nloc), acc_nhat2(1:acc_nloc),         &
    !$acc               acc_nhat3(1:acc_nloc))
    !$acc update device(acc_dia(1:acc_nloc), acc_ct(1:acc_nloc),              &
    !$acc               acc_ud(1:acc_nloc), acc_udt(1:acc_nloc),              &
    !$acc               acc_fn(1:acc_nloc))
    !$acc update device(acc_disk(1:acc_nloc), acc_uc(1:acc_nloc),             &
    !$acc               acc_vc(1:acc_nloc), acc_wc(1:acc_nloc))
end if
if (acc_total_nodes > 0) then
    !$acc update device(acc_i(1:acc_total_nodes), acc_j(1:acc_total_nodes),   &
    !$acc               acc_k(1:acc_total_nodes))
    !$acc update device(acc_ind(1:acc_total_nodes),                           &
    !$acc               acc_ind_t(1:acc_total_nodes),                         &
    !$acc               acc_e1(1:acc_total_nodes), acc_e2(1:acc_total_nodes), &
    !$acc               acc_e3(1:acc_total_nodes))
end if
if (need_alloc) then
    !$acc update device(acc_w_uv, acc_fza_uv)
end if
!$acc update device(acc_use_rotation, acc_inv_tip_speed_ratio)
!$acc update device(acc_nloc, acc_total_nodes)

turbines_acc_ready = .true.

end subroutine turbines_acc_metadata_init

!*******************************************************************************
subroutine turbines_acc_finalize()
!*******************************************************************************
implicit none

if (allocated(acc_offset)) deallocate(acc_offset)
if (allocated(acc_count)) deallocate(acc_count)
if (allocated(acc_i)) deallocate(acc_i)
if (allocated(acc_j)) deallocate(acc_j)
if (allocated(acc_k)) deallocate(acc_k)
if (allocated(acc_icp)) deallocate(acc_icp)
if (allocated(acc_jcp)) deallocate(acc_jcp)
if (allocated(acc_kcp)) deallocate(acc_kcp)
if (allocated(acc_center)) deallocate(acc_center)
if (allocated(acc_ind)) deallocate(acc_ind)
if (allocated(acc_ind_t)) deallocate(acc_ind_t)
if (allocated(acc_e1)) deallocate(acc_e1)
if (allocated(acc_e2)) deallocate(acc_e2)
if (allocated(acc_e3)) deallocate(acc_e3)
if (allocated(acc_nhat1)) deallocate(acc_nhat1)
if (allocated(acc_nhat2)) deallocate(acc_nhat2)
if (allocated(acc_nhat3)) deallocate(acc_nhat3)
if (allocated(acc_dia)) deallocate(acc_dia)
if (allocated(acc_ct)) deallocate(acc_ct)
if (allocated(acc_ud)) deallocate(acc_ud)
if (allocated(acc_udt)) deallocate(acc_udt)
if (allocated(acc_fn)) deallocate(acc_fn)
if (allocated(acc_disk)) deallocate(acc_disk)
if (allocated(acc_uc)) deallocate(acc_uc)
if (allocated(acc_vc)) deallocate(acc_vc)
if (allocated(acc_wc)) deallocate(acc_wc)
if (allocated(acc_w_uv)) deallocate(acc_w_uv)
if (allocated(acc_fza_uv)) deallocate(acc_fza_uv)
acc_nloc = 0
acc_total_nodes = 0
acc_nloc_capacity = 0
acc_node_capacity = 0
turbines_acc_ready = .false.

end subroutine turbines_acc_finalize

!*******************************************************************************
subroutine turbines_acc_sync_device_field(F)
!*******************************************************************************
!
! GPU-aware equivalent of the legacy slab halo sync for a device-resident
! ld-by-ny-by-(lbz:nz) field. It exchanges the same physical planes as
! mpi_sync_real_array(..., lbz, DOWN/DOWNUP), but sends full ld*ny slabs to
! avoid non-contiguous slice temporaries for padded arrays.
!
#ifdef PPMPI
use mpi, only : mpi_sendrecv
#endif
implicit none

real(rprec), dimension(ld,ny,lbz:nz), intent(inout) :: F

#ifdef PPMPI
if (nproc > 1) then
    !$acc host_data use_device(F)
    call mpi_sendrecv(F(1,1,1),  ld*ny, MPI_RPREC, down, 1,                   &
                      F(1,1,nz), ld*ny, MPI_RPREC, up,   1, comm, status, ierr)
    if (lbz == 0) then
        call mpi_sendrecv(F(1,1,nz-1), ld*ny, MPI_RPREC, up,   2,             &
                          F(1,1,0),    ld*ny, MPI_RPREC, down, 2, comm, status, ierr)
    end if
    !$acc end host_data
    if (ierr /= 0) call error(mod_name // '.turbines_acc_sync_device_field',   &
        'MPI slab sync failed with code:', ierr)
end if
#endif

end subroutine turbines_acc_sync_device_field

!*******************************************************************************
subroutine turbines_forcing_acc()
!*******************************************************************************
!
! OpenACC fast path for static, single-rank legacy drag-disk turbines.
!
use sim_param, only : u, v, w, fxa, fya, fza
use functions, only : linear_interp
#ifdef PPMPI
use mpi, only : MPI_Allreduce, MPI_SUM
#endif
implicit none

integer :: node, node0, node1, ii, jj, kk, fid
real(rprec) :: local_disk, local_uc, local_vc, local_wc, ind2
real(rprec) :: top_sum
#ifdef PPMPI
real(rprec), allocatable :: reduce_send(:), reduce_recv(:)
#endif

if (.not. turbines_acc_ready) then
    call turbines_forcing()
    return
end if

if (dyn_theta1 .or. dyn_theta2) then
    do s = 1, nloc
        if (dyn_theta1) wind_farm%turbine(s)%theta1 =                           &
            linear_interp(theta1_time, theta1_arr(s,:), total_time_dim)
        if (dyn_theta2) wind_farm%turbine(s)%theta2 =                           &
            linear_interp(theta2_time, theta2_arr(s,:), total_time_dim)
        if (dyn_Ct_prime) wind_farm%turbine(s)%Ct_prime =                       &
            linear_interp(Ct_prime_time, Ct_prime_arr(s,:), total_time_dim)
    end do
    call turbines_nodes()
    call turbines_acc_metadata_init()
else if (dyn_Ct_prime) then
    do s = 1, acc_nloc
        wind_farm%turbine(s)%Ct_prime =                                       &
            linear_interp(Ct_prime_time, Ct_prime_arr(s,:), total_time_dim)
        acc_ct(s) = wind_farm%turbine(s)%Ct_prime
    end do
    !$acc update device(acc_ct)
end if

! Interpolate w to the uv grid on device.
!$acc parallel loop collapse(3) default(present)
do kk = 1, nz-1
do jj = 1, ny
do ii = 1, ld
    acc_w_uv(ii,jj,kk) = 0.5_rprec * (w(ii,jj,kk+1) + w(ii,jj,kk))
end do
end do
end do

#ifdef PPMPI
if (coord == nproc - 1) then
#endif
    !$acc parallel loop collapse(2) default(present)
    do jj = 1, ny
    do ii = 1, ld
        acc_w_uv(ii,jj,nz) = acc_w_uv(ii,jj,nz-1)
    end do
    end do
#ifdef PPMPI
end if
#endif

call turbines_acc_sync_device_field(acc_w_uv)

! Compute disk-averaged velocity and center velocity on device.
!$acc parallel loop gang default(present) private(node,node0,node1,local_disk,local_uc,local_vc,local_wc,ii,jj,kk)
do s = 1, acc_nloc
    local_disk = 0._rprec
    local_uc = 0._rprec
    local_vc = 0._rprec
    local_wc = 0._rprec
    node0 = acc_offset(s)
    node1 = node0 + acc_count(s) - 1
    !$acc loop seq
    do node = node0, node1
        ii = acc_i(node)
        jj = acc_j(node)
        kk = acc_k(node)
        local_disk = local_disk + dx*dy*dz*acc_ind(node)                      &
            * (acc_nhat1(s)*u(ii,jj,kk) + acc_nhat2(s)*v(ii,jj,kk)            &
            + acc_nhat3(s)*acc_w_uv(ii,jj,kk))
    end do
    if (acc_center(s) /= 0) then
        local_uc = u(acc_icp(s), acc_jcp(s), acc_kcp(s))
        local_vc = v(acc_icp(s), acc_jcp(s), acc_kcp(s))
        local_wc = acc_w_uv(acc_icp(s), acc_jcp(s), acc_kcp(s))
    end if
    acc_disk(s) = local_disk
    acc_uc(s) = local_uc
    acc_vc(s) = local_vc
    acc_wc(s) = local_wc
end do

! Single-rank path: no MPI allreduce is needed.
if (acc_nloc > 0) then
    !$acc update self(acc_disk(1:acc_nloc), acc_uc(1:acc_nloc),               &
    !$acc             acc_vc(1:acc_nloc), acc_wc(1:acc_nloc))
end if

! Reduce small per-turbine scalar diagnostics across ranks. The heavy fields
! stay on device; only O(n_turbines) arrays move through host MPI.
#ifdef PPMPI
if (nproc > 1) then
    allocate(reduce_send(4*acc_nloc), reduce_recv(4*acc_nloc))
    reduce_send(1:acc_nloc) = acc_disk
    reduce_send(acc_nloc+1:2*acc_nloc) = acc_uc
    reduce_send(2*acc_nloc+1:3*acc_nloc) = acc_vc
    reduce_send(3*acc_nloc+1:4*acc_nloc) = acc_wc
    call MPI_Allreduce(reduce_send, reduce_recv, 4*acc_nloc, MPI_rprec,       &
        MPI_SUM, comm, ierr)
    if (ierr /= 0) call error(mod_name // '.turbines_forcing_acc',             &
        'MPI_Allreduce failed with code:', ierr)
    acc_disk(1:acc_nloc) = reduce_recv(1:acc_nloc)
    acc_uc(1:acc_nloc) = reduce_recv(acc_nloc+1:2*acc_nloc)
    acc_vc(1:acc_nloc) = reduce_recv(2*acc_nloc+1:3*acc_nloc)
    acc_wc(1:acc_nloc) = reduce_recv(3*acc_nloc+1:4*acc_nloc)
    deallocate(reduce_send, reduce_recv)
end if
#endif

! Update running disk averages and thrust scalars on host; copy only the small
! per-turbine scalars back to device for force scatter.
if (T_avg_dim > 0.) then
    eps = (dt_dim / T_avg_dim) / (1. + dt_dim / T_avg_dim)
else
    eps = 1.
end if

do s = 1, acc_nloc
    wind_farm%turbine(s)%u_d = acc_disk(s)
    if (adm_correction) then
        wind_farm%turbine(s)%u_d = wind_farm%turbine(s)%u_d                   &
            / (1._rprec + 0.25_rprec                                          &
            * (1._rprec - wind_farm%turbine(s)%turb_ind_func%M) * acc_ct(s))
    end if
    wind_farm%turbine(s)%u_d_T = (1._rprec - eps) * wind_farm%turbine(s)%u_d_T &
        + eps * wind_farm%turbine(s)%u_d
    wind_farm%turbine(s)%f_n = -0.5_rprec * acc_ct(s)                         &
        * abs(wind_farm%turbine(s)%u_d_T) * wind_farm%turbine(s)%u_d_T        &
        * 0.25_rprec * pi * acc_dia(s)**2
    acc_ud(s) = wind_farm%turbine(s)%u_d
    acc_udt(s) = wind_farm%turbine(s)%u_d_T
    acc_fn(s) = wind_farm%turbine(s)%f_n

    if (modulo(jt_total, tbase) == 0 .and. coord == 0) then
        write(forcing_fid(s), *) total_time_dim, acc_uc(s), acc_vc(s),         &
            acc_wc(s), -acc_ud(s), -acc_udt(s), wind_farm%turbine(s)%theta1,   &
            wind_farm%turbine(s)%theta2, acc_ct(s)
    end if
end do

if (acc_nloc > 0) then
    !$acc update device(acc_ud(1:acc_nloc), acc_udt(1:acc_nloc),              &
    !$acc               acc_fn(1:acc_nloc))
end if

! Reset and scatter applied forces on device.
!$acc parallel loop collapse(3) default(present)
do kk = lbz, nz
do jj = 1, ny
do ii = 1, ld
    fxa(ii,jj,kk) = 0._rprec
    fya(ii,jj,kk) = 0._rprec
    fza(ii,jj,kk) = 0._rprec
    acc_fza_uv(ii,jj,kk) = 0._rprec
end do
end do
end do

if (acc_use_rotation /= 0) then
    do s = 1, acc_nloc
        node0 = acc_offset(s)
        node1 = node0 + acc_count(s) - 1
        !$acc parallel loop default(present) private(ii,jj,kk,ind2)
        do node = node0, node1
            ii = acc_i(node)
            jj = acc_j(node)
            kk = acc_k(node)
            ind2 = acc_ind(node)
            fxa(ii,jj,kk) = acc_fn(s) * acc_nhat1(s) * ind2
            fya(ii,jj,kk) = acc_fn(s) * acc_nhat2(s) * ind2
            acc_fza_uv(ii,jj,kk) = acc_fn(s) * acc_nhat3(s) * ind2
            ind2 = acc_ind_t(node) * acc_inv_tip_speed_ratio
            fxa(ii,jj,kk) = fxa(ii,jj,kk) + acc_fn(s) * acc_e1(node) * ind2
            fya(ii,jj,kk) = fya(ii,jj,kk) + acc_fn(s) * acc_e2(node) * ind2
            acc_fza_uv(ii,jj,kk) = acc_fza_uv(ii,jj,kk)                       &
                + acc_fn(s) * acc_e3(node) * ind2
    end do
    end do
else
    do s = 1, acc_nloc
        node0 = acc_offset(s)
        node1 = node0 + acc_count(s) - 1
        !$acc parallel loop default(present) private(ii,jj,kk,ind2)
        do node = node0, node1
            ii = acc_i(node)
            jj = acc_j(node)
            kk = acc_k(node)
            ind2 = acc_ind(node)
            fxa(ii,jj,kk) = acc_fn(s) * acc_nhat1(s) * ind2
            fya(ii,jj,kk) = acc_fn(s) * acc_nhat2(s) * ind2
            acc_fza_uv(ii,jj,kk) = acc_fn(s) * acc_nhat3(s) * ind2
    end do
    end do
end if

call turbines_acc_sync_device_field(fxa)
call turbines_acc_sync_device_field(fya)
call turbines_acc_sync_device_field(acc_fza_uv)

! Interpolate force onto the w grid. The lower overlap/bogus layer remains zero,
! matching the host helper behavior for lbz.
!$acc parallel loop collapse(3) default(present)
do kk = lbz+1, nz
do jj = 1, ny
do ii = 1, ld
    fza(ii,jj,kk) = 0.5_rprec * (acc_fza_uv(ii,jj,kk-1) + acc_fza_uv(ii,jj,kk))
end do
end do
end do

! Spatially average velocity at the top of the domain and write to file.
if (coord .eq. nproc-1) then
    top_sum = 0._rprec
    !$acc parallel loop collapse(2) reduction(+:top_sum) present(u)
    do jj = 1, ny
    do ii = 1, nx
        top_sum = top_sum + u(ii,jj,nz-1)
    end do
    end do
    open(newunit=fid, file=vel_top_dat, status='unknown', form='formatted',    &
        action='write', position='append')
    write(fid,*) total_time, top_sum/(nx*ny)
    close(fid)
end if

end subroutine turbines_forcing_acc
#else
!*******************************************************************************
logical function turbines_acc_available()
!*******************************************************************************
implicit none

turbines_acc_available = .false.

end function turbines_acc_available

!*******************************************************************************
subroutine turbines_forcing_acc()
!*******************************************************************************
implicit none

call turbines_forcing()

end subroutine turbines_forcing_acc
#endif

!*******************************************************************************
subroutine turbines_forcing()
!*******************************************************************************
!
! This subroutine applies the drag-disk forcing
!
use param, only : pi, lbz
use sim_param, only : u, v, w, fxa, fya, fza
use functions, only : linear_interp, interp_to_uv_grid, interp_to_w_grid
use mpi, only : MPI_Allreduce, MPI_SUM
implicit none

character(*), parameter :: sub_name = mod_name // '.turbines_forcing'

real(rprec), pointer :: p_u_d => null(), p_u_d_T => null(), p_f_n => null()
real(rprec), pointer :: p_Ct_prime => null()
integer, pointer :: p_icp => null(), p_jcp => null(), p_kcp => null()

integer :: fid

real(rprec) :: ind2
real(rprec), dimension(nloc) :: disk_avg_vel
real(rprec), dimension(nloc) :: u_vel_center, v_vel_center, w_vel_center
real(rprec), allocatable, dimension(:,:,:) :: w_uv
real(rprec), pointer, dimension(:) :: y, z
real(rprec), dimension(nloc) :: buffer_array

nullify(y,z)
y => grid % y
z => grid % z

allocate(w_uv(ld,ny,lbz:nz))

#ifdef PPMPI
!syncing intermediate w-velocities
call mpi_sync_real_array(w, 0, MPI_SYNC_DOWNUP)
#endif

w_uv = interp_to_uv_grid(w, lbz)

! Do interpolation for dynamically changing parameters
do s = 1, nloc
    if (dyn_theta1) wind_farm%turbine(s)%theta1 =                              &
        linear_interp(theta1_time, theta1_arr(s,:), total_time_dim)
    if (dyn_theta2) wind_farm%turbine(s)%theta2 =                              &
        linear_interp(theta2_time, theta2_arr(s,:), total_time_dim)
    if (dyn_Ct_prime) wind_farm%turbine(s)%Ct_prime =                          &
        linear_interp(Ct_prime_time, Ct_prime_arr(s,:), total_time_dim)
end do

! Recompute the turbine position if theta1 or theta2 can change
if (dyn_theta1 .or. dyn_theta2) call turbines_nodes

!Each processor calculates the weighted disk-averaged velocity
disk_avg_vel = 0._rprec
u_vel_center = 0._rprec
v_vel_center = 0._rprec
w_vel_center = 0._rprec
do s=1,nloc
    ! Calculate total disk-averaged velocity for each turbine
    ! (current, instantaneous) in the normal direction by integrating the
    ! velocity times the indicator function
    do l=1,wind_farm%turbine(s)%num_nodes
        i2 = wind_farm%turbine(s)%nodes(l,1)
        j2 = wind_farm%turbine(s)%nodes(l,2)
        k2 = wind_farm%turbine(s)%nodes(l,3)
        disk_avg_vel(s) = disk_avg_vel(s)                                      &
            + dx*dy*dz*wind_farm%turbine(s)%ind(l)                             &
            * ( wind_farm%turbine(s)%nhat(1)*u(i2,j2,k2)                       &
            + wind_farm%turbine(s)%nhat(2)*v(i2,j2,k2)                         &
            + wind_farm%turbine(s)%nhat(3)*w_uv(i2,j2,k2) )
    end do

    ! Set pointers
    p_icp => wind_farm%turbine(s)%icp
    p_jcp => wind_farm%turbine(s)%jcp
    p_kcp => wind_farm%turbine(s)%kcp

    ! Calculate disk center velocity
    if (wind_farm%turbine(s)%center_in_proc) then
        u_vel_center(s) = u(p_icp, p_jcp, p_kcp)
        v_vel_center(s) = v(p_icp, p_jcp, p_kcp)
        w_vel_center(s) = w_uv(p_icp, p_jcp, p_kcp)
    end if
end do

! Calculate disk velocities by summing all processors and multiplying by disk volume
#ifdef PPMPI
call MPI_Allreduce(disk_avg_vel, buffer_array, nloc, MPI_rprec, MPI_SUM, comm, ierr)
disk_avg_vel = buffer_array
call MPI_Allreduce(u_vel_center, buffer_array, nloc, MPI_rprec, MPI_SUM, comm, ierr)
u_vel_center = buffer_array
call MPI_Allreduce(v_vel_center, buffer_array, nloc, MPI_rprec, MPI_SUM, comm, ierr)
v_vel_center = buffer_array
call MPI_Allreduce(w_vel_center, buffer_array, nloc, MPI_rprec, MPI_SUM, comm, ierr)
w_vel_center = buffer_array
#endif

! Update epsilon for the new timestep (for cfl_dt)
if (T_avg_dim > 0.) then
    eps = (dt_dim / T_avg_dim) / (1. + dt_dim / T_avg_dim)
else
    eps = 1.
end if

! Calculate and apply disk force
do s=1,nloc
    !set pointers
    p_u_d => wind_farm%turbine(s)%u_d
    p_u_d_T => wind_farm%turbine(s)%u_d_T
    p_f_n => wind_farm%turbine(s)%f_n
    p_Ct_prime => wind_farm%turbine(s)%Ct_prime

    !add this current value to the "running average" (first order filter)
    p_u_d = disk_avg_vel(s)
    if (adm_correction) then
        p_u_d = p_u_d /(1 + 0.25_rprec                                         &
            * (1-wind_farm%turbine(s)%turb_ind_func%M)*p_Ct_prime)
    end if
    p_u_d_T = (1.-eps)*p_u_d_T + eps*p_u_d

    !calculate total thrust force for each turbine  (per unit mass)
    !force is normal to the surface (calc from u_d_T, normal to surface)
    !write force to array that will be transferred via MPI
    p_f_n = -0.5*p_Ct_prime*abs(p_u_d_T)*p_u_d_T*0.25*pi*wind_farm%turbine(s)%dia**2

    !write values to file
    if (modulo (jt_total, tbase) == 0 .and. coord == 0) then
        write( forcing_fid(s), *) total_time_dim, u_vel_center(s),         &
            v_vel_center(s), w_vel_center(s), -p_u_d, -p_u_d_T,            &
            wind_farm%turbine(s)%theta1, wind_farm%turbine(s)%theta2,      &
            p_Ct_prime
    end if


    do l=1,wind_farm%turbine(s)%num_nodes
        i2 = wind_farm%turbine(s)%nodes(l,1)
        j2 = wind_farm%turbine(s)%nodes(l,2)
        k2 = wind_farm%turbine(s)%nodes(l,3)
        ind2 = wind_farm%turbine(s)%ind(l)
        fxa(i2,j2,k2) = p_f_n*wind_farm%turbine(s)%nhat(1)*ind2
        fya(i2,j2,k2) = p_f_n*wind_farm%turbine(s)%nhat(2)*ind2
        fza(i2,j2,k2) = p_f_n*wind_farm%turbine(s)%nhat(3)*ind2
        if (use_rotation) then
            ind2 = wind_farm%turbine(s)%ind_t(l)
            fxa(i2,j2,k2) = fxa(i2,j2,k2)                                      &
                + p_f_n*wind_farm%turbine(s)%e_theta(l,1)*ind2/tip_speed_ratio
            fya(i2,j2,k2) = fya(i2,j2,k2)                                      &
                + p_f_n*wind_farm%turbine(s)%e_theta(l,2)*ind2/tip_speed_ratio
            fza(i2,j2,k2) = fza(i2,j2,k2)                                      &
                + p_f_n*wind_farm%turbine(s)%e_theta(l,3)*ind2/tip_speed_ratio
        end if
    end do
end do

! Interpolate force onto the w grid
call mpi_sync_real_array( fxa(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( fya(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( fza(1:nx,1:ny,lbz:nz), 0, MPI_SYNC_DOWNUP )
fza = interp_to_w_grid(fza,lbz)

!spatially average velocity at the top of the domain and write to file
if (coord .eq. nproc-1) then
    open(newunit=fid, file=vel_top_dat, status='unknown', form='formatted',    &
        action='write', position='append')
    write(fid,*) total_time, sum(u(:,:,nz-1))/(nx*ny)
    close(fid)
end if

! Cleanup
deallocate(w_uv)
nullify(y,z)
nullify(p_icp, p_jcp, p_kcp)

end subroutine turbines_forcing

!*******************************************************************************
subroutine turbines_finalize ()
!*******************************************************************************
implicit none

character (*), parameter :: sub_name = mod_name // '.turbines_finalize'

! Persist disk-averaged turbine velocities and averaging time for restart or
! multi-run continuation.
call turbines_checkpoint

#ifdef PPLES_GPU
call turbines_acc_finalize()
#endif

!deallocate
deallocate(wind_farm%turbine)

end subroutine turbines_finalize

!*******************************************************************************
subroutine turbines_checkpoint ()
!*******************************************************************************
!
!
!
implicit none

character (*), parameter :: sub_name = mod_name // '.turbines_checkpoint'
integer :: fid

! Write disk-averaged turbine velocities with T_avg_dim so a continuation run
! can preserve the turbine averaging state.
if (coord == 0) then
    open(newunit=fid, file=u_d_T_dat, status='unknown', form='formatted',      &
        position='rewind')
    do i=1,nloc
        write(fid,*) wind_farm%turbine(i)%u_d_T
    end do
    write(fid,*) T_avg_dim
    close (fid)
end if

end subroutine turbines_checkpoint

!*******************************************************************************
subroutine turbine_vel_init(zo_high)
!*******************************************************************************
!
! called from ic.f90 if initu, lbc_mom==1, S_FLAG are all false.
! this accounts for the turbines when creating the initial velocity profile.
!
use param, only: zo
implicit none
character (*), parameter :: sub_name = mod_name // '.turbine_vel_init'

real(rprec), intent(inout) :: zo_high
real(rprec) :: cft, nu_w, exp_KE, induction_factor, Ct_noprime

! Convert Ct' to Ct
! a = Ct'/(4+Ct'), Ct = 4a(1-a)
induction_factor = Ct_prime / (4._rprec + Ct_prime)
Ct_noprime = 4*(induction_factor) * (1 - induction_factor)

! friction coefficient, cft
cft = pi*Ct_noprime/(4.*sx*sy)

!wake viscosity
nu_w = 28.*sqrt(0.5*cft)

!turbine friction height, Calaf, Phys. Fluids 22, 2010
zo_high = height_all*(1.+0.5*dia_all/height_all)**(nu_w/(1.+nu_w))* &
  exp(-1.*(0.5*cft/(vonk**2) + (log(height_all/zo* &
  (1.-0.5*dia_all/height_all)**(nu_w/(1.+nu_w))) )**(-2) )**(-0.5) )

exp_KE =  0.5*(log(0.45/zo_high)/0.4)**2

if(.false.) then
    write(*,*) 'sx,sy,cft: ',sx,sy,cft
    write(*,*) 'nu_w: ',nu_w
    write(*,*) 'zo_high: ',zo_high
    write(*,*) 'approx expected KE: ', exp_KE
end if
end subroutine turbine_vel_init

!*******************************************************************************
subroutine place_turbines
!*******************************************************************************
!
! This subroutine places the turbines on the domain. It also sets the values for
! each individual turbine. After the subroutine is called, the following values
! are set for each turbine in wind_farm: xloc, yloc, height, dia, thk, theta1,
! theta2, and Ct_prime.
!
use param, only: pi, z_i
use messages, only : error, warn
implicit none

character(*), parameter :: sub_name = mod_name // '.place_turbines'

real(rprec) :: sxx, syy, shift_base, const
real(rprec) :: dummy, dummy2
logical :: exst
integer :: fid

! Read parameters from file if needed
if (read_param) then
    ! Check if file exists and open
    inquire (file = param_dat, exist = exst)
    if (.not. exst) then
        call error (sub_name, 'file ' // param_dat // 'does not exist')
    end if

    ! Check that there are enough lines from which to read data
    nloc = count_lines(param_dat)
    if (nloc < num_x*num_y) then
        nloc = num_x*num_y
        call error(sub_name, param_dat // 'must have num_x*num_y lines')
    else if (nloc > num_x*num_y) then
        call warn(sub_name, param_dat // ' has more than num_x*num_y lines. '  &
                  // 'Only reading first num_x*num_y lines')
    end if

    ! Read from parameters file, which should be in this format:
    ! xloc [meters], yloc [meters], height [meters], dia [meters], thk [meters],
    ! theta1 [degrees], theta2 [degrees], Ct_prime [-]
    write(*,*) "Reading from", param_dat
    open(newunit=fid, file=param_dat, status='unknown', form='formatted',      &
        position='rewind')
    do k = 1, nloc
        read(fid,*) wind_farm%turbine(k)%xloc, wind_farm%turbine(k)%yloc,      &
            wind_farm%turbine(k)%height, wind_farm%turbine(k)%dia,             &
            wind_farm%turbine(k)%thk, wind_farm%turbine(k)%theta1,             &
            wind_farm%turbine(k)%theta2, wind_farm%turbine(k)%Ct_prime
    end do
    close(fid)

    ! Make lengths dimensionless
    do k = 1, nloc
        wind_farm%turbine(k)%xloc = wind_farm%turbine(k)%xloc / z_i
        wind_farm%turbine(k)%yloc = wind_farm%turbine(k)%yloc / z_i
        wind_farm%turbine(k)%height = wind_farm%turbine(k)%height / z_i
        wind_farm%turbine(k)%dia = wind_farm%turbine(k)%dia / z_i
        wind_farm%turbine(k)%thk = wind_farm%turbine(k)%thk / z_i
    end do
else
    ! Set values for each turbine based on values in input file
    wind_farm%turbine(:)%height = height_all
    wind_farm%turbine(:)%dia = dia_all
    wind_farm%turbine(:)%thk = thk_all
    wind_farm%turbine(:)%theta1 = theta1_all
    wind_farm%turbine(:)%theta2 = theta2_all
    wind_farm%turbine(:)%Ct_prime = Ct_prime

    ! Set baseline locations (evenly spaced, not staggered aka aligned)
    k = 1
    sxx = sx * dia_all  ! x-spacing with units to match those of L_x
    syy = sy * dia_all  ! y-spacing
    do i = 1,num_x
        do j = 1,num_y
            wind_farm%turbine(k)%xloc = sxx*real(2*i-1)/2
            wind_farm%turbine(k)%yloc = syy*real(2*j-1)/2
            k = k + 1
        end do
    end do

    ! Place turbines based on orientation flag
    ! This will shift the placement relative to the baseline locations abive
    select case (orientation)
        ! Evenly-spaced, not staggered
        case (1)

        ! Evenly-spaced, horizontally staggered only
        ! Shift each row according to stag_perc
        case (2)
            do i = 2, num_x
                do k = 1+num_y*(i-1), num_y*i
                    shift_base = syy * stag_perc/100.
                    wind_farm%turbine(k)%yloc = mod( wind_farm%turbine(k)%yloc &
                                                    + (i-1)*shift_base , L_y )
                end do
            end do

        ! Evenly-spaced, only vertically staggered (by rows)
        case (3)
            ! Make even rows taller
            do i = 2, num_x, 2
                do k = 1+num_y*(i-1), num_y*i
                    wind_farm%turbine(k)%height = height_all*(1.+stag_perc/100.)
                end do
            end do
            ! Make odd rows shorter
            do i = 1, num_x, 2
                do k = 1+num_y*(i-1), num_y*i
                    wind_farm%turbine(k)%height = height_all*(1.-stag_perc/100.)
                end do
            end do

        ! Evenly-spaced, only vertically staggered, checkerboard pattern
        case (4)
            k = 1
            do i = 1, num_x
                do j = 1, num_y
                    ! this should alternate between 1, -1
                    const = 2.*mod(real(i+j),2.)-1.
                    wind_farm%turbine(k)%height = height_all                   &
                                                  *(1.+const*stag_perc/100.)
                    k = k + 1
                end do
            end do

        ! Aligned, but shifted forward for efficient use of simulation space
        ! during CPS runs
        case (5)
        ! Shift in spanwise direction: Note that stag_perc is now used
            k=1
            dummy=stag_perc                                                    &
                  *(wind_farm%turbine(2)%yloc - wind_farm%turbine(1)%yloc)
            do i = 1, num_x
                do j = 1, num_y
                    dummy2=dummy*(i-1)
                    wind_farm%turbine(k)%yloc=mod( wind_farm%turbine(k)%yloc   &
                                                  + dummy2,L_y)
                    k=k+1
                end do
            end do

        case default
            call error (sub_name, 'invalid orientation')

    end select
end if

end subroutine place_turbines

!*******************************************************************************
subroutine read_control_files
!*******************************************************************************
!
! This subroutine reads the input files for dynamic controls with theta1,
! theta2, and Ct_prime. This is calles from turbines_init.
!
use param, only: pi
implicit none

character(*), parameter :: sub_name = mod_name // '.place_turbines'

integer :: fid, i, num_t

! Read the theta1 input data
if (dyn_theta1) then
    ! Count number of entries and allocate
    num_t = count_lines(theta1_dat)
    allocate( theta1_time(num_t) )
    allocate( theta1_arr(nloc, num_t) )

    ! Read values from file
    open(newunit=fid, file=theta1_dat, status='unknown', form='formatted',     &
        position='rewind')
    do i = 1, num_t
        read(fid,*) theta1_time(i), theta1_arr(:,i)
    end do
end if

! Read the theta2 input data
if (dyn_theta2) then
    ! Count number of entries and allocate
    num_t = count_lines(theta2_dat)
    allocate( theta2_time(num_t) )
    allocate( theta2_arr(nloc, num_t) )

    ! Read values from file
    open(newunit=fid, file=theta2_dat, status='unknown', form='formatted',     &
        position='rewind')
    do i = 1, num_t
        read(fid,*) theta2_time(i), theta2_arr(:,i)
    end do
end if

! Read the Ct_prime input data
if (dyn_Ct_prime) then
    ! Count number of entries and allocate
    num_t = count_lines(Ct_prime_dat)
    allocate( Ct_prime_time(num_t) )
    allocate( Ct_prime_arr(nloc, num_t) )

    ! Read values from file
    open(newunit=fid, file=Ct_prime_dat, status='unknown', form='formatted',   &
        position='rewind')
    do i = 1, num_t
        read(fid,*) Ct_prime_time(i), Ct_prime_arr(:,i)
    end do
end if

end subroutine read_control_files

end module turbines
