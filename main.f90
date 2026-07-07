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
!!

!*******************************************************************************
program main
!*******************************************************************************
!
! Main file for lesgo solver
! Contains main time loop
!
! Navigation map:
!   - setup: input parsing, initialization, optional turbine/scalar setup
!   - main timestep loop: derivatives, wall stress, SGS, convection, forcing,
!     pressure, projection, diagnostics, and output
!   - GPU/CPU boundary points: explicit device updates around wall, SGS, ATM,
!     and pressure stages
!   - shutdown: final output, checkpoints, module finalizers
!   - local helpers: main_read_env_real and main_cuda_sync near file end
!
! Keep this file as orchestration code.  Module-specific kernels and policy
! logic belong in their owning modules unless the dependency is truly global.

use types, only : rprec
use clock_m, only : clock_t
use param
use sim_param
use io, only : energy, output_loop, output_final, jt_total
use io, only : write_tau_wall_bot, write_tau_wall_top
use derivatives, only : filt_da_vel, ddz_vel
use cfl_util, only : get_cfl_dt, get_max_cfl
use sgs_stag_util, only : sgs_stag
use forcing, only : forcing_random, forcing_applied, forcing_induced, project
use functions, only: get_tau_wall_bot, get_tau_wall_top

use iwmles, only : iwm_lhs_update_due
use sgs_param, only : SGS_MODEL_SMAGORINSKY, SGS_MODEL_STANDARD_DYNAMIC,      &
                      SGS_MODEL_SCALE_DEP_DYNAMIC,                            &
                      SGS_MODEL_LAGRANGE_SIMILARITY,                          &
                      SGS_MODEL_LAGRANGE_SCALE_DEP

#ifdef PPMPI
use mpi
use mpi_defs, only : mpi_sync_real_array, MPI_SYNC_DOWN
use cuda_mpi_debug, only : mpi_dbg_sendrecv_r
#endif

#ifdef PPLVLSET
use level_set, only : level_set_global_CA, level_set_vel_err
use level_set_base, only : global_CA_calc
#if defined(PPSGS_GPU) && defined(PPLES_GPU)
use sgs_param, only : S11, S12, S13, S22, S23, S33, Nu_t, Cs_opt2,            &
                      F_LM, F_MM, F_QN, F_NN, Beta, Tn_all
#ifdef PPDYN_TN
use sgs_param, only : F_ee2, F_deedt2, ee_past
#endif
#endif
#endif

#if defined(PPATM) && defined(PPLES_GPU)
use atm_lesgo_interface, only : atm_lesgo_forcing
#endif

#ifdef PPTURBINES
use turbines, only : turbines_forcing, turbine_vel_init
#endif

#ifdef PPSCALARS
use scalars, only : buoyancy_force, scalars_transport, scalars_deriv,        &
    passive_scalar
#endif

use sponge, only : sponge_force
use coriolis, only : coriolis_calc, coriolis_forcing, alpha, G, phi_actual
use messages, only : error, mesg

#ifdef PPCONVEC_GPU
use convec_gpu_m, only : convec_gpu
#endif

#ifdef PPDERIVS_GPU
use derivatives_gpu_m, only : filt_da_gpu, ddz_uv_gpu, ddz_w_gpu
#endif

#ifdef PPPRESS_GPU
use press_gpu_m, only : press_stag_array_gpu
#endif

#ifdef PPSGS_GPU
use sgs_gpu_m, only : sgs_stag_gpu, divstress_uv_gpu, divstress_w_gpu
#endif


implicit none

character (*), parameter :: prog_name = 'main'
integer :: nca
character(:), allocatable :: ca

integer :: jt_step, nstart
integer :: jx, jy, jz
real(rprec) :: rmsdivvel, ke, maxcfl, tt

type(clock_t) :: clock, clock_total, clock_forcing

! --- --------------------------------------------------------------
type(clock_t) :: clock_derivs, clock_sgs, clock_convec
type(clock_t) :: clock_press, clock_project, clock_output
type(clock_t) :: clock_derivs_xy, clock_derivs_z
type(clock_t) :: clock_sgs_model, clock_sgs_halo
type(clock_t) :: clock_sgs_divuv, clock_sgs_divw

! ---
real(rprec) :: total_time_derivs = 0.0_rprec
real(rprec) :: total_time_derivs_xy = 0.0_rprec
real(rprec) :: total_time_derivs_z = 0.0_rprec
real(rprec) :: total_time_sgs = 0.0_rprec
real(rprec) :: total_time_sgs_model = 0.0_rprec
real(rprec) :: total_time_sgs_halo = 0.0_rprec
real(rprec) :: total_time_sgs_divuv = 0.0_rprec
real(rprec) :: total_time_sgs_divw = 0.0_rprec
real(rprec) :: total_time_convec = 0.0_rprec
real(rprec) :: total_time_press = 0.0_rprec
real(rprec) :: total_time_project = 0.0_rprec
real(rprec) :: total_time_output = 0.0_rprec
real(rprec) :: total_time_named = 0.0_rprec
real(rprec) :: total_time_other = 0.0_rprec
real(rprec) :: cpu_ref_time_total_runtime = 0.0_rprec
real(rprec) :: cpu_ref_time_forcing_runtime = 0.0_rprec
real(rprec) :: cpu_ref_time_other_runtime = 0.0_rprec
logical :: cpu_ref_time_total_available = .false.

real(rprec), parameter :: cpu_ref_time_derivs = 0.6389535E-01_rprec
real(rprec), parameter :: cpu_ref_time_sgs = 0.2967275E+00_rprec
real(rprec), parameter :: cpu_ref_time_convec = 0.2444816E+00_rprec
real(rprec), parameter :: cpu_ref_time_press = 0.1355035E+00_rprec
real(rprec), parameter :: cpu_ref_time_project = 0.5635220E-01_rprec
real(rprec), parameter :: cpu_ref_time_named = cpu_ref_time_derivs +       &
    cpu_ref_time_sgs + cpu_ref_time_convec + cpu_ref_time_press +          &
    cpu_ref_time_project

! Measure total time in forcing function
real(rprec) :: clock_total_f = 0.0

#ifdef PPMPI
! Buffers used for MPI communication
real(rprec) :: rbuffer
real(rprec) :: maxdummy ! Used to calculate maximum with mpi_allreduce
real(rprec) :: tau_top   ! Used to write top wall stress at first proc
#endif

! Initialize MPI
#ifdef PPMPI
call mpi_init (ierr)
#endif

! Get path if needed
nca = COMMAND_ARGUMENT_COUNT()
if (nca == 1) then
    call GET_COMMAND_ARGUMENT(1, length=nca)
    allocate(character(nca) :: ca)
    call GET_COMMAND_ARGUMENT(1, value=ca)
    allocate(path, source = './' // ca // '/')
else
    allocate(path, source='./')
endif

! Start the clocks, both local and total
call clock%start

! Initialize time variable
tt = 0
jt = 0
jt_total = 0

! Initialize all data
call initialize()

if(coord == 0) then
    call clock%stop
#ifdef PPMPI
    write(*,'(1a,E15.7)') 'Initialization wall time: ', clock % time
#else
    write(*,'(1a,E15.7)') 'Initialization cpu time: ', clock % time
#endif
endif

call clock_total%start

! Initialize starting loop index
! If new simulation jt_total=0 by definition, if restarting jt_total
! provided by total_time.dat
nstart = jt_total+1

! BEGIN TIME LOOP
time_loop: do jt_step = nstart, nsteps

    ! Get the starting time for the iteration
    call clock%start

    if (use_cfl_dt) then

        dt_f = dt
        dt = get_cfl_dt()
        dt_dim = dt * z_i / u_star

        tadv1 = 1._rprec + 0.5_rprec * dt / dt_f
        tadv2 = 1._rprec - tadv1

    end if

   ! Advance time
   jt_total = jt_step
   jt = jt + 1
   total_time = total_time + dt
   total_time_dim = total_time_dim + dt_dim
   tt = tt+dt

    ! Save previous time's right-hand-sides for Adams-Bashforth Integration
    ! NOTE: RHS does not contain the pressure gradient
#ifdef PPLES_GPU
    !$acc parallel loop collapse(3) default(present) async(1)
    do jz = lbz, nz
    do jy = 1, ny
    do jx = 1, ld
        RHSx_f(jx,jy,jz) = RHSx(jx,jy,jz)
        RHSy_f(jx,jy,jz) = RHSy(jx,jy,jz)
        RHSz_f(jx,jy,jz) = RHSz(jx,jy,jz)
    end do
    end do
    end do
#else
    RHSx_f = RHSx
    RHSy_f = RHSy
    RHSz_f = RHSz
#endif

    ! ------------------------------------------------------
    !/// CALCULATE DERIVATIVES                          ///
    call clock_derivs%start

    ! Calculate velocity derivatives
    ! Calculate dudx, dudy, dvdx, dvdy, dwdx, dwdy (in Fourier space)
    call clock_derivs_xy%start
#ifdef PPDERIVS_GPU
    call filt_da_gpu(u, dudx, dudy, lbz)
    call filt_da_gpu(v, dvdx, dvdy, lbz)
    call filt_da_gpu(w, dwdx, dwdy, lbz)
#else
    call filt_da_vel(u, v, w, dudx, dudy, dvdx, dvdy, dwdx, dwdy, lbz)
#endif
    call clock_derivs_xy%stop
    total_time_derivs_xy = clock_derivs_xy%time

    ! Calculate dudz, dvdz using finite differences (for 1:nz on uv-nodes)
    !  except bottom coord, only 2:nz
    call clock_derivs_z%start
#ifdef PPDERIVS_GPU
    call ddz_uv_gpu(u, dudz, lbz)
    call ddz_uv_gpu(v, dvdz, lbz)
    call ddz_w_gpu(w, dwdz, lbz)
    ! No active PPSGS_GPU SGS model uses the old CPU SGS fallback here.
    ! Wallstress only reads u,v at the wall plane(s), so sync only those.
    if (.false.) then
        !$acc wait(1)
        !$acc update self(u, v, w, dudx, dudy, dudz, dvdx, dvdy, dvdz,          &
        !$acc             dwdx, dwdy, dwdz)
    else
#if defined(PPSGS_GPU)
#ifdef PPSCALARS
        ! In scalar builds, lbc_mom=2 equilibrium and lbc_mom=3 IWM are device-resident.
        ! Keep host staging for scalar DNS/free/upper wallstress paths.
        if (coord == 0 .and. lbc_mom /= 3 .and.                                &
            lbc_mom /= 2) then
            !$acc wait(1)
            !$acc update self(u(:,:,1), v(:,:,1))
        end if
        if (coord == nproc-1) then
            !$acc wait(1)
            !$acc update self(u(:,:,nz-1), v(:,:,nz-1))
        end if
#else
        ! wallstress runs on the device in this path, including the OpenACC
        ! integral wall model for lbc_mom=3, so no wall-plane host staging is
        ! needed here.
#endif
#else
        if (coord == 0) then
            !$acc wait(1)
            if (lbc_mom == 3) then
                if (iwm_lhs_update_due()) then
                    !$acc update self(u(:,:,1), v(:,:,1), w(:,:,2), p(:,:,1))
                end if
            else
                !$acc update self(u(:,:,1), v(:,:,1))
            end if
        end if
        if (coord == nproc-1) then
            !$acc wait(1)
            !$acc update self(u(:,:,nz-1), v(:,:,nz-1))
        end if
#endif
    end if
#else
    call ddz_vel(u, v, w, dudz, dvdz, dwdz, lbz)
#endif
    call clock_derivs_z%stop
    total_time_derivs_z = clock_derivs_z%time

    call clock_derivs%stop
    total_time_derivs =  clock_derivs%time

#ifdef PPSCALARS
    call scalars_deriv()
#endif

    ! Calculate wall stress and derivatives at the wall
    ! (txz, tyz, dudz, dvdz at jz=1)
    ! using the velocity log-law
    ! MPI: bottom and top processes only
    if (coord == 0 .or. coord == nproc-1) then
        call wallstress()
#ifdef PPLES_GPU
#if defined(PPSGS_GPU)
#ifdef PPSCALARS
        ! For scalar builds, lower equilibrium and IWM wallstress write these
        ! planes on device. Other scalar wallstress paths still use host.
        if (coord == 0 .and. lbc_mom /= 3 .and.                                &
            lbc_mom /= 2) then
            !$acc update device(dudz(:,:,1), dvdz(:,:,1), txz(:,:,1), tyz(:,:,1))
        end if
        if (coord == nproc-1) then
            !$acc update device(dudz(:,:,nz), dvdz(:,:,nz), txz(:,:,nz), tyz(:,:,nz))
        end if
#else
        ! Device-resident wallstress writes txz/tyz/dudz/dvdz wall planes on
        ! the GPU directly, including IWM lbc_mom=3.
#endif
#else
        if (coord == 0) then
            !$acc update device(dudz(:,:,1), dvdz(:,:,1), txz(:,:,1), tyz(:,:,1))
        end if
        if (coord == nproc-1) then
            !$acc update device(dudz(:,:,nz), dvdz(:,:,nz), txz(:,:,nz), tyz(:,:,nz))
        end if
#endif
#endif
    end if

    ! --- Calculate subgrid stress

    call clock_sgs%start

    ! Calculate turbulent (subgrid) stress for entire domain
    !   using the model specified in param.f90 (Smag, LASD, etc)
    !   MPI: txx, txy, tyy, tzz at 1:nz-1; txz, tyz at 1:nz (stress-free lid)
    call clock_sgs_model%start
#ifdef PPSGS_GPU
    if (sgs .and. sgs_model /= SGS_MODEL_SMAGORINSKY .and.                   &
        sgs_model /= SGS_MODEL_STANDARD_DYNAMIC .and.                         &
        sgs_model /= SGS_MODEL_SCALE_DEP_DYNAMIC .and.                        &
        sgs_model /= SGS_MODEL_LAGRANGE_SIMILARITY .and.                     &
        sgs_model /= SGS_MODEL_LAGRANGE_SCALE_DEP) then
        call error(prog_name,                                                   &
            'PPSGS_GPU currently supports active SGS models 1, 2, 3, 4, and 5 only')
    end if
#if defined(PPLVLSET) && defined(PPLES_GPU)
    ! LVLSET has no OpenACC SGS path yet. Use the established CPU LVLSET SGS
    ! routines and explicitly bridge the resident LES/SGS arrays.
    !$acc wait(1)
    !$acc update self(u, v, w, dudx, dudy, dudz, dvdx, dvdy, dvdz,             &
    !$acc             dwdx, dwdy, dwdz, txx, txy, txz, tyy, tyz, tzz,         &
    !$acc             S11, S12, S13, S22, S23, S33, Nu_t, Cs_opt2,            &
    !$acc             F_LM, F_MM, F_QN, F_NN, Beta, Tn_all)
#ifdef PPDYN_TN
    !$acc update self(F_ee2, F_deedt2, ee_past)
#endif
    call sgs_stag()
    !$acc update device(txx, txy, txz, tyy, tyz, tzz,                         &
    !$acc               S11, S12, S13, S22, S23, S33, Nu_t, Cs_opt2,          &
    !$acc               F_LM, F_MM, F_QN, F_NN, Beta, Tn_all)
#ifdef PPDYN_TN
    !$acc update device(F_ee2, F_deedt2, ee_past)
#endif
#else
    call sgs_stag_gpu()
#endif
#else
    call sgs_stag()
#endif
    call clock_sgs_model%stop
    total_time_sgs_model = clock_sgs_model%time


    ! Exchange ghost node information (since coords overlap) for tau_zz
    !   send info up (from nz-1 below to 0 above)
    total_time_sgs_halo = 0.0_rprec
#ifdef PPMPI
#ifdef PPSGS_GPU
    call clock_sgs_halo%start
    !$acc wait(1)
#ifdef PPGPU_AWARE_MPI
    !$acc host_data use_device(tzz)
    call mpi_sendrecv (tzz(1,1,nz-1), ld*ny, MPI_RPREC, up, 6,                &
                       tzz(1,1,0), ld*ny, MPI_RPREC, down, 6,                 &
                       comm, status, ierr)
    !$acc end host_data
#else
    !$acc update self(tzz(:,:,nz-1))
    call mpi_sendrecv (tzz(:,:,nz-1), ld*ny, MPI_RPREC, up, 6,                &
                       tzz(:,:,0), ld*ny, MPI_RPREC, down, 6,                 &
                       comm, status, ierr)
    !$acc update device(tzz(:,:,0))
#endif
    call clock_sgs_halo%stop
    total_time_sgs_halo = clock_sgs_halo%time
#else
    call clock_sgs_halo%start
    call mpi_dbg_sendrecv_r (tzz(1,1,nz-1), ld*ny, MPI_RPREC, up, 6,          &
                       tzz(1,1,0), ld*ny, MPI_RPREC, down, 6,                 &
                       comm, status, ierr, 'main_tzz_halo')
    call clock_sgs_halo%stop
    total_time_sgs_halo = clock_sgs_halo%time
#endif
#endif

    ! Compute divergence of SGS shear stresses
    ! the divt's and the diagonal elements of t are not equivalenced
    ! in this version. Provides divtz 1:nz-1, except 1:nz at top process
    call clock_sgs_divuv%start
#ifdef PPSGS_GPU
#if defined(PPLVLSET) && defined(PPLES_GPU)
    ! LVLSET stress boundary treatment is still host-only. Keep div-stress on
    ! the same path as the LVLSET SGS bridge, then refresh device divtau.
    call divstress_uv(divtx, divty, txx, txy, txz, tyy, tyz)
    !$acc update device(divtx, divty)
#else
    call divstress_uv_gpu(divtx, divty)
#endif
#else
    call divstress_uv(divtx, divty, txx, txy, txz, tyy, tyz)
#endif
    call clock_sgs_divuv%stop
    total_time_sgs_divuv = clock_sgs_divuv%time

    call clock_sgs_divw%start
#ifdef PPSGS_GPU
#if defined(PPLVLSET) && defined(PPLES_GPU)
    call divstress_w(divtz, txz, tyz, tzz)
    !$acc update device(divtz)
#else
    call divstress_w_gpu(divtz)
#endif
#else
    call divstress_w(divtz, txz, tyz, tzz)
#endif
    call clock_sgs_divw%stop
    total_time_sgs_divw = clock_sgs_divw%time

    ! ----------
    call clock_sgs%stop
    total_time_sgs =  clock_sgs%time

    ! --- Calculate convection term
    call clock_convec%start

    ! Calculates u x (omega) term in physical space. Uses 3/2 rule for
    ! dealiasing. Stores this term in RHS (right hand side) variable
#ifdef PPCONVEC_GPU
    call convec_gpu()
#else
    call convec()
#endif
    ! -------------------------
    call clock_convec%stop
    total_time_convec = clock_convec%time

#if defined(PPATM) && defined(PPLES_GPU)
    ! ATM phase 1 (w_uv interp, blade update, device velocity sampling, host
    ! blade-force model, MPI gather) runs HERE - convec_gpu now returns with
    ! its async(1) kernels still draining (~60 ms backlog), so phase 1's
    ! ~30 ms of host work overlaps it. Same inputs (u, v, w are not modified
    ! until the velocity update) and the same internal order of ATM
    ! operations - just earlier within the step. Phase 2 (convolution into
    ! the device fxa/fya/fza, Cl correction, apply, output) stays at the
    ! original forcing_applied call site below.
    call clock_forcing%start
    call atm_lesgo_forcing(phase=1)
    call clock_forcing%stop
    clock_total_f = clock_total_f + clock_forcing % time
#endif

    ! Add div-tau term to RHS variable
    !   this will be used for pressure calculation
#ifdef PPLES_GPU
    !$acc parallel loop collapse(3) default(present) async(1)
    do jz = 1, nz - 1
    do jy = 1, ny
    do jx = 1, ld
        RHSx(jx,jy,jz) = -RHSx(jx,jy,jz) - divtx(jx,jy,jz)
        RHSy(jx,jy,jz) = -RHSy(jx,jy,jz) - divty(jx,jy,jz)
        RHSz(jx,jy,jz) = -RHSz(jx,jy,jz) - divtz(jx,jy,jz)
    end do
    end do
    end do
    if (coord == nproc-1) then
        !$acc parallel loop collapse(2) default(present) async(1)
        do jy = 1, ny
        do jx = 1, ld
            RHSz(jx,jy,nz) = -RHSz(jx,jy,nz) - divtz(jx,jy,nz)
        end do
        end do
    end if
#else
    RHSx(:,:,1:nz-1) = -RHSx(:,:,1:nz-1) - divtx(:,:,1:nz-1)
    RHSy(:,:,1:nz-1) = -RHSy(:,:,1:nz-1) - divty(:,:,1:nz-1)
    RHSz(:,:,1:nz-1) = -RHSz(:,:,1:nz-1) - divtz(:,:,1:nz-1)
    if (coord == nproc-1) RHSz(:,:,nz) = -RHSz(:,:,nz)-divtz(:,:,nz)
#endif

    call coriolis_calc()

#ifdef PPSCALARS
    call scalars_transport()
    call buoyancy_force()
#endif
    call sponge_force()

    !--calculate u^(*) (intermediate vel field)
    !  at this stage, p, dpdx_i are from previous time step
    !  (assumes old dpdx has NOT been added to RHSx_f, etc)
    !  we add force (mean press forcing) here so that u^(*) is as close
    !  to the final velocity as possible
    if (use_mean_p_force) then
#ifdef PPLES_GPU
        !$acc parallel loop collapse(3) default(present) async(1)
        do jz = 1, nz - 1
        do jy = 1, ny
        do jx = 1, ld
            RHSx(jx,jy,jz) = RHSx(jx,jy,jz) + mean_p_force_x
            RHSy(jx,jy,jz) = RHSy(jx,jy,jz) + mean_p_force_y
        end do
        end do
        end do
#else
        RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + mean_p_force_x
        RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) + mean_p_force_y
#endif
    end if

    ! Optional random forcing, i.e. to help prevent relaminarization
    if (use_random_force .and. jt_total < stop_random_force) then
        call forcing_random()
    end if

    !//////////////////////////////////////////////////////
    !/// APPLIED FORCES                                 ///
    !//////////////////////////////////////////////////////
    !  In order to save memory the arrays fxa, fya, and fza are now only defined when needed.
    !  For Levelset RNS all three arrays are assigned.
    !  For turbines at the moment only fxa is assigned.
    !  Look in forcing_applied for calculation of forces.
    !  Look in sim_param.f90 for the assignment of the arrays.

    !  Applied forcing (forces are added to RHS{x,y,z})

    ! Calculate forcing time
    call clock_forcing%start

    ! Apply forcing. These forces will later go into RHS
    call forcing_applied()

    ! Calculate forcing time
    call clock_forcing%stop

    ! Calculate the total time of the forcing
    clock_total_f = clock_total_f + clock_forcing % time

    !  Update RHS with applied forcing
#if defined(PPTURBINES) || defined(PPATM)
#ifdef PPLES_GPU
    !$acc parallel loop collapse(3) default(present) async(1)
    do jz = 1, nz - 1
    do jy = 1, ny
    do jx = 1, ld
        RHSx(jx,jy,jz) = RHSx(jx,jy,jz) + fxa(jx,jy,jz)
        RHSy(jx,jy,jz) = RHSy(jx,jy,jz) + fya(jx,jy,jz)
        RHSz(jx,jy,jz) = RHSz(jx,jy,jz) + fza(jx,jy,jz)
    end do
    end do
    end do
#else
    RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + fxa(:,:,1:nz-1)
    RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) + fya(:,:,1:nz-1)
    RHSz(:,:,1:nz-1) = RHSz(:,:,1:nz-1) + fza(:,:,1:nz-1)
#endif
#endif

    !//////////////////////////////////////////////////////
    !/// EULER INTEGRATION CHECK                        ///
    !//////////////////////////////////////////////////////
    ! Set RHS*_f if necessary (first timestep)
    if ((jt_total == 1) .and. (.not. initu)) then
        ! if initu, then this is read from the initialization file
        ! else for the first step put RHS_f=RHS
        !--i.e. at first step, take an Euler step
#ifdef PPLES_GPU
        !$acc parallel loop collapse(3) default(present) async(1)
        do jz = lbz, nz
        do jy = 1, ny
        do jx = 1, ld
            RHSx_f(jx,jy,jz) = RHSx(jx,jy,jz)
            RHSy_f(jx,jy,jz) = RHSy(jx,jy,jz)
            RHSz_f(jx,jy,jz) = RHSz(jx,jy,jz)
        end do
        end do
        end do
#else
        RHSx_f = RHSx
        RHSy_f = RHSy
        RHSz_f = RHSz
#endif
    end if

    !//////////////////////////////////////////////////////
    !/// INTERMEDIATE VELOCITY                          ///
    !//////////////////////////////////////////////////////
    ! Calculate intermediate velocity field
    !   only 1:nz-1 are valid
#ifdef PPLES_GPU
    !$acc parallel loop collapse(3) default(present) async(1)
    do jz = 1, nz - 1
    do jy = 1, ny
    do jx = 1, ld
        u(jx,jy,jz) = u(jx,jy,jz) + dt *                                  &
            (tadv1 * RHSx(jx,jy,jz) + tadv2 * RHSx_f(jx,jy,jz))
        v(jx,jy,jz) = v(jx,jy,jz) + dt *                                  &
            (tadv1 * RHSy(jx,jy,jz) + tadv2 * RHSy_f(jx,jy,jz))
        w(jx,jy,jz) = w(jx,jy,jz) + dt *                                  &
            (tadv1 * RHSz(jx,jy,jz) + tadv2 * RHSz_f(jx,jy,jz))
    end do
    end do
    end do
    if (coord == nproc-1) then
        !$acc parallel loop collapse(2) default(present) async(1)
        do jy = 1, ny
        do jx = 1, ld
            w(jx,jy,nz) = w(jx,jy,nz) + dt *                              &
                (tadv1 * RHSz(jx,jy,nz) + tadv2 * RHSz_f(jx,jy,nz))
        end do
        end do
    end if
#else
    u(:,:,1:nz-1) = u(:,:,1:nz-1) +                                            &
        dt * ( tadv1 * RHSx(:,:,1:nz-1) + tadv2 * RHSx_f(:,:,1:nz-1) )
    v(:,:,1:nz-1) = v(:,:,1:nz-1) +                                            &
        dt * ( tadv1 * RHSy(:,:,1:nz-1) + tadv2 * RHSy_f(:,:,1:nz-1) )
    w(:,:,1:nz-1) = w(:,:,1:nz-1) +                                            &
        dt * ( tadv1 * RHSz(:,:,1:nz-1) + tadv2 * RHSz_f(:,:,1:nz-1) )
    if (coord == nproc-1) then
        w(:,:,nz) = w(:,:,nz) +                                                &
            dt * ( tadv1 * RHSz(:,:,nz) + tadv2 * RHSz_f(:,:,nz) )
    end if
#endif

    ! Set unused values to BOGUS so unintended uses will be noticable
#ifdef PPSAFETYMODE
#ifdef PPLES_GPU
#ifdef PPMPI
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
    do jx = 1, ld
        u(jx,jy,0) = BOGUS
        v(jx,jy,0) = BOGUS
        w(jx,jy,0) = BOGUS
    end do
    end do
#endif
    !$acc parallel loop collapse(2) default(present) async(1)
    do jy = 1, ny
    do jx = 1, ld
        u(jx,jy,nz) = BOGUS
        v(jx,jy,nz) = BOGUS
        if(coord < nproc-1) w(jx,jy,nz) = BOGUS
    end do
    end do
#else
#ifdef PPMPI
    u(:,:,0) = BOGUS
    v(:,:,0) = BOGUS
    w(:,:,0) = BOGUS
#endif
    u(:,:,nz) = BOGUS
    v(:,:,nz) = BOGUS
    if(coord < nproc-1) w(:,:,nz) = BOGUS
#endif
#endif

    !//////////////////////////////////////////////////////
    !/// PRESSURE SOLUTION                              ///
    !//////////////////////////////////////////////////////
    ! Solve Poisson equation for pressure
    !   div of momentum eqn + continuity (div-vel=0) yields Poisson eqn
    !   do not need to store p --> only need gradient
    !   provides p, dpdx, dpdy at 0:nz-1 and dpdz at 1:nz-1

    ! --- Calculate pressure and pressure gradients
    call clock_press%start

#ifdef PPPRESS_GPU
#if defined(PPLVLSET) && defined(PPLES_GPU)
    ! LVLSET bridge: pressure still has to consume host RHS values produced by
    ! the host LVLSET/divstress fallback, then return pressure gradients to the
    ! device for the remaining explicit-residency update/projection kernels.
    !$acc wait(1)
    !$acc update self(u, v, w, divtz, RHSx, RHSy, RHSz)
    call press_stag_array()
    !$acc update device(RHSx, RHSy, RHSz, p, dpdx, dpdy, dpdz)
#else
    call press_stag_array_gpu()
#endif
#else
    call press_stag_array()
#endif

    ! Add pressure gradients to RHS variables (for next time step)
    !   could avoid storing pressure gradients - add directly to RHS
#ifdef PPLES_GPU
    !$acc parallel loop collapse(3) default(present) async(1)
    do jz = 1, nz - 1
    do jy = 1, ny
    do jx = 1, ld
        RHSx(jx,jy,jz) = RHSx(jx,jy,jz) - dpdx(jx,jy,jz)
        RHSy(jx,jy,jz) = RHSy(jx,jy,jz) - dpdy(jx,jy,jz)
        RHSz(jx,jy,jz) = RHSz(jx,jy,jz) - dpdz(jx,jy,jz)
    end do
    end do
    end do
    if(coord == nproc-1) then
        !$acc parallel loop collapse(2) default(present) async(1)
        do jy = 1, ny
        do jx = 1, ld
            RHSz(jx,jy,nz) = RHSz(jx,jy,nz) - dpdz(jx,jy,nz)
        end do
        end do
    end if
#else
    RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) - dpdx(:,:,1:nz-1)
    RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) - dpdy(:,:,1:nz-1)
    RHSz(:,:,1:nz-1) = RHSz(:,:,1:nz-1) - dpdz(:,:,1:nz-1)
    if(coord == nproc-1) then
        RHSz(:,:,nz) = RHSz(:,:,nz) - dpdz(:,:,nz)
    end if
#endif
    ! -------------------------------
    call clock_press%stop
    total_time_press =  clock_press%time

    !//////////////////////////////////////////////////////
    !/// INDUCED FORCES                                 ///
    !//////////////////////////////////////////////////////
    ! Calculate external forces induced forces. These are
    ! stored in fx,fy,fz arrays. We are calling induced
    ! forces before applied forces as some of the applied
    ! forces (RNS) depend on the induced forces and the
    ! two are assumed independent
    call forcing_induced()

    !//////////////////////////////////////////////////////
    !/// PROJECTION STEP                                ///
    !//////////////////////////////////////////////////////
    ! Projection method provides u,v,w for jz=1:nz
    !   uses fx,fy,fz calculated above
    !   for MPI: syncs 1 -> Nz and Nz-1 -> 0 nodes info for u,v,w

    ! --------------------
    call clock_project%start

    call project ()
    ! --------------------
    call clock_project%stop
    total_time_project = clock_project%time

    ! Write ke to file
    if (modulo (jt_total, nenergy) == 0) then
#if defined(PPLES_GPU)
        ! energy() reduces on the device now - the 765 MB u/v/w D2H this
        ! used to do every nenergy steps is gone.
        !$acc wait(1)
#elif defined(PPLES_GPU)
        !$acc wait(1)
        !$acc update self(u, v, w)
#endif
        call energy(ke)
    end if

#ifdef PPLVLSET
#if defined(PPLES_GPU)
    if (global_CA_calc) then
        !$acc wait(1)
        !$acc update self(u, fx, fy, fz)
    end if
#endif
    if (global_CA_calc) call level_set_global_CA()
#endif

    ! Write output files
    call output_loop()

    ! Check the total time of the simulation up to this point on the master
    ! node and send this to all

    if (modulo (jt_total, wbase) == 0) then

        ! Get the ending time for the iteration
        call clock%stop
        call clock_total%stop

        ! Calculate rms divergence of velocity
        ! only written to screen, not used otherwise
#if defined(PPDERIVS_GPU) && defined(PPLES_GPU)
        ! rmsdiv() reduces on the device now - no dudx/dvdy/dwdz D2H.
        !$acc wait(1)
#elif defined(PPDERIVS_GPU)
        !$acc wait(1)
        !$acc update self(dudx, dvdy, dwdz)
#endif
#if defined(PPSGS_GPU)
        ! get_tau_wall_bot/top reduce txz/tyz directly on the OpenACC device,
        ! so wbase diagnostics no longer need wall-plane device-to-host copies.
        !$acc wait(1)
#endif
        call rmsdiv(rmsdivvel)
        maxcfl = get_max_cfl()

        ! This takes care of the clock times, to obtain the quantities based
        ! on all the processors, not just processor 0
#ifdef PPMPI
        call mpi_allreduce(clock % time, maxdummy, 1, mpi_rprec,               &
            MPI_MAX, comm, ierr)
        clock % time = maxdummy
        call mpi_allreduce(clock_total % time, maxdummy, 1, mpi_rprec,         &
            MPI_MAX, comm, ierr)
        clock_total % time = maxdummy
        call mpi_allreduce(clock_forcing % time, maxdummy, 1, mpi_rprec,       &
            MPI_MAX, comm, ierr)
        clock_forcing % time = maxdummy
        call mpi_allreduce(clock_total_f , maxdummy, 1, mpi_rprec,             &
            MPI_MAX, comm, ierr)
        clock_total_f = maxdummy
#endif

        ! Send top wall stress to bottom process
#ifdef PPMPI
        if (coord == nproc-1) then
            tau_top = get_tau_wall_top()
        else
            tau_top = 0._rprec
        endif

        call mpi_allreduce(tau_top, maxdummy, 1, mpi_rprec,               &
            MPI_SUM, comm, ierr)
        tau_top = maxdummy
#endif

        total_time_named = total_time_derivs + total_time_sgs +             &
            total_time_convec + total_time_press + total_time_project
        total_time_other = max(clock % time - total_time_named -            &
            clock_forcing % time, 0.0_rprec)

        cpu_ref_time_total_runtime = 0.0_rprec
        cpu_ref_time_forcing_runtime = 0.0_rprec
        cpu_ref_time_other_runtime = -1.0_rprec
        cpu_ref_time_total_available = .false.

        call main_read_env_real('LESGO_CPU_REF_TIME_TOTAL',                 &
            cpu_ref_time_total_runtime, cpu_ref_time_total_available)
        call main_read_env_real('LESGO_CPU_REF_TIME_FORCING',               &
            cpu_ref_time_forcing_runtime)

        if (cpu_ref_time_total_available) then
            cpu_ref_time_other_runtime = max(cpu_ref_time_total_runtime -    &
                cpu_ref_time_named - cpu_ref_time_forcing_runtime,           &
                0.0_rprec)
        endif

            if (coord == 0) then
            write(*,*)
            write(*,'(a)') '==================================================='
            write(*,'(a)') 'Time step information:'
            write(*,'(a,i9)') '  Iteration: ', jt_total
            write(*,'(a,E15.7)') '  Time step: ', dt
            write(*,'(a,E15.7)') '  Dimensional time: ', total_time_dim
            write(*,'(a,E15.7)') '  CFL: ', maxcfl
            write(*,'(a,2E15.7)') '  AB2 TADV1, TADV2: ', tadv1, tadv2
            write(*,*)
            write(*,'(a)') 'Flow field information:'
            write(*,'(a,E15.7)') '  Velocity divergence metric: ', rmsdivvel
            write(*,'(a,E15.7)') '  Kinetic energy: ', ke
            write(*,'(a,E15.7)') '  Bot wall stress: ', get_tau_wall_bot()
#ifdef PPMPI
            write(*,'(a,E15.7)') '  Top wall stress: ', tau_top
#else
            write(*,'(a,E15.7)') '  Top wall stress: ', get_tau_wall_top()
#endif
            write(*,*)
            write(*,'(1a)') 'Simulation wall times (s): '
            write(*,'(1a,E15.7)') '  Iteration: ', clock % time
            write(*,'(1a,E15.7)') '  Cumulative: ', clock_total % time
            write(*,'(1a,E15.7)') '  Forcing: ', clock_forcing % time
            write(*,'(1a,E15.7)') '  Cumulative Forcing: ', clock_total_f
            write(*,'(1a,E15.7)') '  Forcing %: ',                             &
                clock_total_f /clock_total % time


            write(*,'(a)') '---------------------------------------------------'
            write(*,'(a)') 'Sub-component Cumulative Times (s):'
            write(*,'(1a,E15.7)') '  Derivatives: ', total_time_derivs
            write(*,'(1a,E15.7)') '    Derivatives xy/filter: ',               &
                total_time_derivs_xy
            write(*,'(1a,E15.7)') '    Derivatives z: ', total_time_derivs_z
            write(*,'(1a,E15.7)') '  SGS & Stresses: ', total_time_sgs
            write(*,'(1a,E15.7)') '    SGS model/stress build: ',              &
                total_time_sgs_model
            write(*,'(1a,E15.7)') '    SGS tzz halo: ', total_time_sgs_halo
            write(*,'(1a,E15.7)') '    SGS divstress_uv: ',                    &
                total_time_sgs_divuv
            write(*,'(1a,E15.7)') '    SGS divstress_w: ', total_time_sgs_divw
            write(*,'(1a,E15.7)') '  Convection: ', total_time_convec
            write(*,'(1a,E15.7)') '  Pressure Solver: ', total_time_press
            write(*,'(1a,E15.7)') '  Projection: ', total_time_project
            write(*,'(1a,E15.7)') '  Other: ', total_time_other

            write(*,'(a)') 'CPU Reference vs GPU Time (s, 24 MPI CPU baseline):'
            write(*,'(1a,E15.7,1a,E15.7,1a,E12.4)')                           &
                '  Derivatives: CPU=', cpu_ref_time_derivs, ' GPU=',           &
                total_time_derivs, ' CPU/GPU=',                                &
                cpu_ref_time_derivs/max(total_time_derivs, tiny(1._rprec))
            write(*,'(1a,E15.7,1a,E15.7,1a,E12.4)')                           &
                '  SGS & Stresses: CPU=', cpu_ref_time_sgs, ' GPU=',           &
                total_time_sgs, ' CPU/GPU=',                                   &
                cpu_ref_time_sgs/max(total_time_sgs, tiny(1._rprec))
            write(*,'(1a,E15.7,1a,E15.7,1a,E12.4)')                           &
                '  Convection: CPU=', cpu_ref_time_convec, ' GPU=',            &
                total_time_convec, ' CPU/GPU=',                                &
                cpu_ref_time_convec/max(total_time_convec, tiny(1._rprec))
            write(*,'(1a,E15.7,1a,E15.7,1a,E12.4)')                           &
                '  Pressure Solver: CPU=', cpu_ref_time_press, ' GPU=',        &
                total_time_press, ' CPU/GPU=',                                 &
                cpu_ref_time_press/max(total_time_press, tiny(1._rprec))
            write(*,'(1a,E15.7,1a,E15.7,1a,E12.4)')                           &
                '  Projection: CPU=', cpu_ref_time_project, ' GPU=',           &
                total_time_project, ' CPU/GPU=',                               &
                cpu_ref_time_project/max(total_time_project, tiny(1._rprec))
            if (cpu_ref_time_other_runtime >= 0.0_rprec) then
                write(*,'(1a,E15.7,1a,E15.7,1a,E12.4)')                       &
                    '  Other: CPU=', cpu_ref_time_other_runtime, ' GPU=',      &
                    total_time_other, ' CPU/GPU=',                             &
                    cpu_ref_time_other_runtime/max(total_time_other,           &
                    tiny(1._rprec))
            else
                write(*,'(1a,E15.7,1a)') '  Other: CPU=N/A GPU=',             &
                    total_time_other, ' CPU/GPU=N/A'
            endif



            if (coriolis_forcing > 0) then
                write(*,*)
                write(*,'(1a)') 'Coriolis parameters: '
                write(*,'(1a,E15.7)') '  G: ', G
                write(*,'(1a,E15.7)') '  alpha: ', alpha
                if (coriolis_forcing == 2) then
                    write(*,'(1a,E15.7)') '  wind direction: ', phi_actual
                end if
            end if
            write(*,'(a)') '==================================================='
            call write_tau_wall_bot()
        end if
        if(coord == nproc-1) then
            call write_tau_wall_top()
        end if

#ifdef PPMPI
        call mpi_barrier(comm, ierr)
#endif

        ! Check if we are to check the allowable runtime
        if (runtime > 0) then

#ifdef PPMPI
            ! Determine the processor that has used most time and communicate
            ! this. Needed to make sure that all processors stop at the same
            ! time and not just some of them
            call mpi_allreduce(clock_total % time, rbuffer, 1, MPI_RPREC,      &
                MPI_MAX, MPI_COMM_WORLD, ierr)
            clock_total % time = rbuffer
#endif

            ! If maximum time is surpassed go to the end of the program
            if ( clock_total % time >= real(runtime,rprec) ) then
                call mesg( prog_name,                                          &
                    'Specified runtime exceeded. Exiting simulation.')
                exit time_loop
            endif

       endif

    end if

end do time_loop
! END TIME LOOP

! Finalize
close(2)

! Write total_time.dat and tavg files
call output_final()

! Stop wall clock
call clock_total%stop
#ifdef PPMPI
if (coord == 0) write(*,"(a,e15.7)") 'Simulation wall time (s) : ',            &
    clock_total % time
#else
if (coord == 0) write(*,"(a,e15.7)") 'Simulation cpu time (s) : ',             &
    clock_total % time
#endif

call finalize()

if(coord == 0 ) write(*,'(a)') 'Simulation complete'

contains

!*******************************************************************************
subroutine main_read_env_real(name, value, available)
!*******************************************************************************
implicit none

character(len=*), intent(in) :: name
real(rprec), intent(out) :: value
logical, intent(out), optional :: available
character(len=64) :: env_value
integer :: env_len, env_stat, env_iostat

value = 0.0_rprec
if (present(available)) available = .false.

call get_environment_variable(name, env_value, env_len, env_stat)
if (env_stat == 0 .and. env_len > 0) then
    read(env_value(1:env_len), *, iostat=env_iostat) value
    if (env_iostat /= 0) then
        value = 0.0_rprec
    else if (present(available)) then
        available = .true.
    endif
endif

end subroutine main_read_env_real


end program main
