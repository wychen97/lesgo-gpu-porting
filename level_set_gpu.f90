!!
!! GPU kernels for the immersed-boundary Level Set treatment.
!!
#if defined(PPLVLSET) && defined(PPLVLSET_GPU)
#define LS_GRID_ARGS ld,nx,ny,nz,lbz,dx,dy,dz,L_x,L_y
#define LS_HALO_GRID_ARGS ld,nx,ny,nz,dx,dy,dz,L_x,L_y
module level_set_gpu_m
! Navigation map:
!   1. Persistent workspace and interpolation startup validation
!   2. Stress interpolation, extrapolation, and immersed-boundary treatment
!   3. Desired-velocity and Level Set forcing kernels
!   4. Two-dimensional and three-dimensional smoothing kernels
!   5. Scale-dependent Lagrangian SGS model 4/5 kernels
!   6. Packed MPI halo helpers for geometry, velocity, stress, and SGS state
use types, only : rprec
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
use, intrinsic :: iso_fortran_env, only : int64
use param, only : ld, nx, ny, nz, lbz, dx, dy, dz, L_x, L_y, coord, nproc, &
                  BOGUS
use messages, only : error
use level_set_base, only : phi, physBC, use_smooth_tau, smooth_mode,           &
    use_log_profile, use_extrap_tau_log, use_extrap_tau_simple,               &
    use_modify_dutdn, zo_level_set, nphitop, nphibot, nveltop, nvelbot,       &
    ntautop, ntaubot,                                                          &
    nFMMtop, nFMMbot, phitop, phibot, utop, vtop, wtop, ubot, vbot, wbot,    &
    txxtop, txytop, txztop, tyytop, tyztop, tzztop,                          &
    txxbot, txybot, txzbot, tyybot, tyzbot, tzzbot, FMMtop, FMMbot
use sim_param, only : u, v, w, txx, txy, txz, tyy, tyz, tzz,                  &
    dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
implicit none
private
real(rprec), allocatable, save :: tau_source(:,:,:)
real(rprec), allocatable, save :: fmm_source(:,:,:)
integer, save :: workspace_ld=0, workspace_ny=0, workspace_nz=0
integer, save :: workspace_lbz=0
integer, save :: interp_oob_count=0
logical, save :: interp_bounds_check_enabled=.false.
logical, save :: interp_selftest_active=.false.
!$acc declare create(interp_oob_count,interp_bounds_check_enabled,            &
!$acc                interp_selftest_active)

public :: level_set_bc_gpu_core, level_set_gpu_interp_selftest,             &
          level_set_gpu_workspace_init,                                     &
          level_set_gpu_workspace_finalize,                                 &
          level_set_gpu_workspace_bytes,                                    &
          level_set_gpu_interp_bounds_configure,                             &
          level_set_gpu_interp_bounds_begin,                                 &
          level_set_gpu_interp_bounds_end,                                   &
          interp_tau_gpu, extrap_tau_simple_gpu, smooth_field_gpu,           &
          smooth_tau_gpu, level_set_lag_dyn_gpu_core,                        &
          level_set_desired_velocity_gpu

contains

!*******************************************************************************
subroutine level_set_gpu_interp_bounds_configure()
!*******************************************************************************
character(32) :: setting
integer :: stat

interp_bounds_check_enabled=.false.
call get_environment_variable('LESGO_LVLSET_INTERP_BOUNDS_CHECK',setting,  &
                              status=stat)
if (stat == 0) then
  select case(trim(adjustl(setting)))
    case('1','true','TRUE','True','on','ON','On','yes','YES','Yes')
      interp_bounds_check_enabled=.true.
  end select
end if
interp_oob_count=0
interp_selftest_active=.false.
!$acc update device(interp_oob_count,interp_bounds_check_enabled,             &
!$acc               interp_selftest_active)
end subroutine level_set_gpu_interp_bounds_configure

!*******************************************************************************
subroutine level_set_gpu_interp_bounds_begin()
!*******************************************************************************
if (.not. interp_bounds_check_enabled) return
interp_oob_count=0
!$acc update device(interp_oob_count)
end subroutine level_set_gpu_interp_bounds_begin

!*******************************************************************************
subroutine level_set_gpu_interp_bounds_end(context)
!*******************************************************************************
character(*), intent(in) :: context
if (.not. interp_bounds_check_enabled) return
!$acc wait(1)
!$acc update self(interp_oob_count)
if (interp_oob_count > 0) then
  call error('level_set_gpu_m.' // trim(context),                         &
      'GPU interpolation requested data outside the local/halo contract',&
      interp_oob_count)
end if
end subroutine level_set_gpu_interp_bounds_end

!*******************************************************************************
subroutine ls_record_interp_oob()
!$acc routine seq
!*******************************************************************************
if (interp_bounds_check_enabled) then
  !$acc atomic update
  interp_oob_count=interp_oob_count+1
end if
end subroutine ls_record_interp_oob

!*******************************************************************************
subroutine level_set_gpu_workspace_init(need_tau_source,need_fmm_source)
!*******************************************************************************
! Immutable source snapshots are allocated only for active algorithms.
logical, intent(in) :: need_tau_source,need_fmm_source

if (allocated(tau_source) .or. allocated(fmm_source)) then
  if (workspace_ld == ld .and. workspace_ny == ny .and.                  &
      workspace_nz == nz .and. workspace_lbz == lbz) then
    if (need_tau_source .and. .not. allocated(tau_source)) then
      allocate(tau_source(ld,ny,lbz:nz))
      !$acc enter data create(tau_source)
    else if (.not. need_tau_source .and. allocated(tau_source)) then
      !$acc exit data delete(tau_source)
      deallocate(tau_source)
    end if
    if (need_fmm_source .and. .not. allocated(fmm_source)) then
      allocate(fmm_source(ld,ny,1:nz))
      !$acc enter data create(fmm_source)
    else if (.not. need_fmm_source .and. allocated(fmm_source)) then
      !$acc exit data delete(fmm_source)
      deallocate(fmm_source)
    end if
    return
  else
    call level_set_gpu_workspace_finalize()
  end if
end if
if (need_tau_source) then
  allocate(tau_source(ld,ny,lbz:nz))
  !$acc enter data create(tau_source)
end if
if (need_fmm_source) then
  allocate(fmm_source(ld,ny,1:nz))
  !$acc enter data create(fmm_source)
end if
workspace_ld=ld
workspace_ny=ny
workspace_nz=nz
workspace_lbz=lbz
end subroutine level_set_gpu_workspace_init

!*******************************************************************************
subroutine level_set_gpu_workspace_finalize()
!*******************************************************************************
if (allocated(tau_source)) then
  !$acc exit data delete(tau_source)
  deallocate(tau_source)
end if
if (allocated(fmm_source)) then
  !$acc exit data delete(fmm_source)
  deallocate(fmm_source)
end if
workspace_ld=0
workspace_ny=0
workspace_nz=0
workspace_lbz=0
end subroutine level_set_gpu_workspace_finalize

!*******************************************************************************
integer(int64) function level_set_gpu_workspace_bytes() result(bytes)
!*******************************************************************************
bytes=0_int64
if (allocated(tau_source)) bytes=bytes+size(tau_source,kind=int64)*         &
    int(storage_size(0._rprec)/8,kind=int64)
if (allocated(fmm_source)) bytes=bytes+size(fmm_source,kind=int64)*         &
    int(storage_size(0._rprec)/8,kind=int64)
end function level_set_gpu_workspace_bytes

!*******************************************************************************
subroutine level_set_gpu_interp_selftest(max_error, invalid_count, failed)
!*******************************************************************************
! Exercise device interpolation independently of the production geometry.
! The cases cover affine off-grid points, periodic x/y seams, both z halos,
! and field lower bounds of zero and one.
real(rprec), intent(out) :: max_error
integer, intent(out) :: invalid_count
logical, intent(out) :: failed
integer, parameter :: tnx=4, tny=3, tnz=4, tld=tnx+1, nhalo=2
integer, parameter :: nvalid=10, nsample=11
real(rprec), parameter :: tdx=0.75_rprec, tdy=1.25_rprec, tdz=0.5_rprec
real(rprec), parameter :: tlx=real(tnx,rprec)*tdx
real(rprec), parameter :: tly=real(tny,rprec)*tdy
real(rprec), parameter :: c0=1.25_rprec, cx=0.5_rprec
real(rprec), parameter :: cy=-0.75_rprec, cz=1.1_rprec
real(rprec), allocatable :: a0(:,:,:), a1(:,:,:)
real(rprec), allocatable :: bot0(:,:,:), top0(:,:,:)
real(rprec), allocatable :: bot1(:,:,:), top1(:,:,:)
real(rprec) :: observed(nsample), expected(nvalid), reference_scale
logical :: bounds_guard_failed
integer :: i, j, k, kb, n

allocate(a0(tld,tny,0:tnz), a1(tld,tny,1:tnz))
allocate(bot0(tld,tny,nhalo), top0(tld,tny,nhalo))
allocate(bot1(tld,tny,nhalo), top1(tld,tny,nhalo))

do k=0,tnz
  do j=1,tny
    do i=1,tld
      a0(i,j,k)=affine_value(real(i-1,rprec)*tdx,                         &
          real(j-1,rprec)*tdy,(real(k,rprec)-0.5_rprec)*tdz)
    end do
  end do
end do
do k=1,tnz
  do j=1,tny
    do i=1,tld
      a1(i,j,k)=affine_value(real(i-1,rprec)*tdx,                         &
          real(j-1,rprec)*tdy,real(k-1,rprec)*tdz)
    end do
  end do
end do
do kb=1,nhalo
  do j=1,tny
    do i=1,tld
      k=kb-nhalo-1
      bot0(i,j,kb)=affine_value(real(i-1,rprec)*tdx,                      &
          real(j-1,rprec)*tdy,(real(k,rprec)-0.5_rprec)*tdz)
      k=tnz+kb
      top0(i,j,kb)=affine_value(real(i-1,rprec)*tdx,                      &
          real(j-1,rprec)*tdy,(real(k,rprec)-0.5_rprec)*tdz)
      k=kb-nhalo
      bot1(i,j,kb)=affine_value(real(i-1,rprec)*tdx,                      &
          real(j-1,rprec)*tdy,real(k-1,rprec)*tdz)
      k=tnz+kb
      top1(i,j,kb)=affine_value(real(i-1,rprec)*tdx,                      &
          real(j-1,rprec)*tdy,real(k-1,rprec)*tdz)
    end do
  end do
end do

expected(1)=affine_value(0.37_rprec*tdx,0.42_rprec*tdy,1.15_rprec*tdz)
expected(2)=expected(1)
expected(3)=affine_value(1.2_rprec*tdx,0.33_rprec*tdy,1.4_rprec*tdz)
expected(4)=expected(3)
expected(5)=affine_value(0.8_rprec*tdx,0.6_rprec*tdy,-0.75_rprec*tdz)
expected(6)=affine_value(1.3_rprec*tdx,0.2_rprec*tdy,                       &
                         (real(tnz,rprec)-0.25_rprec)*tdz)
expected(7)=affine_value(0.6_rprec*tdx,0.7_rprec*tdy,-0.25_rprec*tdz)
expected(8)=affine_value(1.1_rprec*tdx,0.4_rprec*tdy,                       &
                         (real(tnz,rprec)-0.5_rprec)*tdz)
expected(9)=affine_value(0.75_rprec*tdx,0.4_rprec*tdy,1.2_rprec*tdz)
expected(10)=expected(9)

interp_selftest_active=.true.
!$acc update device(interp_selftest_active)
!$acc data copyin(a0,a1,bot0,top0,bot1,top1) copyout(observed)
!$acc serial default(present)
observed(1)=ls_interp_field(a0,0.37_rprec*tdx,0.42_rprec*tdy,              &
    1.15_rprec*tdz,.false.,tld,tnx,tny,tnz,0,tdx,tdy,tdz,tlx,tly)
observed(2)=ls_interp_field(a0,0.37_rprec*tdx+tlx,0.42_rprec*tdy-2._rprec*tly,&
    1.15_rprec*tdz,.false.,tld,tnx,tny,tnz,0,tdx,tdy,tdz,tlx,tly)
observed(3)=ls_interp_field(a1,1.2_rprec*tdx,0.33_rprec*tdy,1.4_rprec*tdz, &
    .true.,tld,tnx,tny,tnz,1,tdx,tdy,tdz,tlx,tly)
observed(4)=ls_interp_field(a1,1.2_rprec*tdx-2._rprec*tlx,                &
    0.33_rprec*tdy+tly,1.4_rprec*tdz,.true.,tld,tnx,tny,tnz,1,           &
    tdx,tdy,tdz,tlx,tly)
observed(5)=ls_interp_field_halo(a0,bot0,top0,nhalo,nhalo,0,              &
    0.8_rprec*tdx,0.6_rprec*tdy,-0.75_rprec*tdz,.false.,tld,tnx,tny,tnz,&
    tdx,tdy,tdz,tlx,tly)
observed(6)=ls_interp_field_halo(a0,bot0,top0,nhalo,nhalo,0,              &
    1.3_rprec*tdx,0.2_rprec*tdy,(real(tnz,rprec)-0.25_rprec)*tdz,.false.,&
    tld,tnx,tny,tnz,tdx,tdy,tdz,tlx,tly)
observed(7)=ls_interp_field_halo(a1,bot1,top1,nhalo,nhalo,1,              &
    0.6_rprec*tdx,0.7_rprec*tdy,-0.25_rprec*tdz,.true.,tld,tnx,tny,tnz, &
    tdx,tdy,tdz,tlx,tly)
observed(8)=ls_interp_field_halo(a1,bot1,top1,nhalo,nhalo,1,              &
    1.1_rprec*tdx,0.4_rprec*tdy,(real(tnz,rprec)-0.5_rprec)*tdz,.true.,  &
    tld,tnx,tny,tnz,tdx,tdy,tdz,tlx,tly)
observed(9)=ls_interp_field(a0,-0.25_rprec*tdx,-0.2_rprec*tdy,            &
    1.2_rprec*tdz,.false.,tld,tnx,tny,tnz,0,tdx,tdy,tdz,tlx,tly)
observed(10)=ls_interp_field(a0,tlx-0.25_rprec*tdx,                       &
    tly-0.2_rprec*tdy,1.2_rprec*tdz,.false.,tld,tnx,tny,tnz,0,           &
    tdx,tdy,tdz,tlx,tly)
! Deliberately provide one lower halo plane for a request that needs two.
! The guard must return BOGUS rather than reading before bot0.
observed(11)=ls_interp_field_halo(a0,bot0,top0,1,nhalo,0,                 &
    0.8_rprec*tdx,0.6_rprec*tdy,-1.75_rprec*tdz,.false.,tld,tnx,tny,tnz,&
    tdx,tdy,tdz,tlx,tly)
!$acc end serial
!$acc end data
interp_selftest_active=.false.
interp_oob_count=0
!$acc update device(interp_selftest_active,interp_oob_count)

max_error=0._rprec
invalid_count=0
do n=1,nvalid
  if (.not. ieee_is_finite(observed(n)) .or.                              &
      abs(observed(n)) >= 0.5_rprec*abs(BOGUS)) then
    invalid_count=invalid_count+1
  else
    max_error=max(max_error,abs(observed(n)-expected(n)))
  end if
end do
bounds_guard_failed=.not. ieee_is_finite(observed(nsample)) .or.          &
    abs(observed(nsample)) < 0.5_rprec*abs(BOGUS)
reference_scale=max(1._rprec,maxval(abs(expected)))
failed=invalid_count > 0 .or.                                               &
       max_error > 5000._rprec*epsilon(1._rprec)*reference_scale .or.       &
       bounds_guard_failed

deallocate(a0,a1,bot0,top0,bot1,top1)

contains

pure real(rprec) function affine_value(x,y,z) result(value)
real(rprec), intent(in) :: x,y,z
value=c0+cx*x+cy*y+cz*z
end function affine_value
end subroutine level_set_gpu_interp_selftest

!*******************************************************************************
real(rprec) function ls_interp_field(a, x, y, z, w_node, ldim, nxv, nyv,     &
    nzv, lbzv, dxv, dyv, dzv, lxv, lyv) result(value)
!$acc routine seq
!*******************************************************************************
! Trilinear interpolation on either the uv/phi grid or the staggered w grid.
! Periodicity is arithmetic so device code does not depend on host pointers.
real(rprec), intent(in) :: a(*)
real(rprec), intent(in) :: x, y, z
logical, intent(in) :: w_node
integer, intent(in) :: ldim, nxv, nyv, nzv, lbzv
real(rprec), intent(in) :: dxv, dyv, dzv, lxv, lyv
integer :: i0, i1, j0, j1, k0, k1
integer :: q000, q100, q010, q110, q001, q101, q011, q111
real(rprec) :: xm, ym, ax, ay, az, zs
real(rprec) :: w00, w10, w01, w11

xm = modulo(x, lxv)
ym = modulo(y, lyv)
i0 = min(nxv, max(1, int(floor(xm / dxv)) + 1))
j0 = min(nyv, max(1, int(floor(ym / dyv)) + 1))
i1 = modulo(i0, nxv) + 1
j1 = modulo(j0, nyv) + 1
ax = modulo(xm, dxv) / dxv
ay = modulo(ym, dyv) / dyv

zs = z / dzv
if (w_node) then
  k0 = int(floor(zs + 1._rprec))
  az = modulo(z, dzv) / dzv
else
  k0 = int(floor(zs + 0.5_rprec))
  az = zs - (floor(zs + 0.5_rprec) - 0.5_rprec)
end if
k1 = k0 + 1

! Non-halo interpolation is used at physical boundaries and in non-MPI runs.
! Match the CPU contract instead of flattening an out-of-range array index.
if (k0 < max(1,lbzv) .or. k1 > nzv) then
  call ls_record_interp_oob()
  value=BOGUS
  return
end if

w00 = (1._rprec - ax) * (1._rprec - ay)
w10 = ax * (1._rprec - ay)
w01 = (1._rprec - ax) * ay
w11 = ax * ay
q000=i0+(j0-1)*ldim+(k0-lbzv)*ldim*nyv
q100=i1+(j0-1)*ldim+(k0-lbzv)*ldim*nyv
q010=i0+(j1-1)*ldim+(k0-lbzv)*ldim*nyv
q110=i1+(j1-1)*ldim+(k0-lbzv)*ldim*nyv
q001=i0+(j0-1)*ldim+(k1-lbzv)*ldim*nyv
q101=i1+(j0-1)*ldim+(k1-lbzv)*ldim*nyv
q011=i0+(j1-1)*ldim+(k1-lbzv)*ldim*nyv
q111=i1+(j1-1)*ldim+(k1-lbzv)*ldim*nyv
value=(1._rprec-az)*(w00*a(q000)+w10*a(q100)+w01*a(q010)+w11*a(q110)) +    &
      az*(w00*a(q001)+w10*a(q101)+w01*a(q011)+w11*a(q111))
end function ls_interp_field

!*******************************************************************************
real(rprec) function ls_halo_value(a,abot,atop,nbot,ntop,albzv,k,i,j,       &
    ldim,nyv,nzv) result(value)
!$acc routine seq
!*******************************************************************************
real(rprec), intent(in) :: a(*), abot(*), atop(*)
integer, intent(in) :: nbot,ntop,albzv,k,i,j,ldim,nyv,nzv
integer :: kb, q

if (.not. interp_selftest_active) then
  if ((coord == 0 .and. k < 1) .or.                                  &
      (coord == nproc-1 .and. k > nzv)) then
    call ls_record_interp_oob()
    value=BOGUS
    return
  end if
end if

if (k < albzv) then
  kb=nbot+k+1-albzv
  if (kb < 1 .or. kb > nbot) then
    call ls_record_interp_oob()
    value=BOGUS
    return
  end if
  q=i+(j-1)*ldim+(kb-1)*ldim*nyv
  value=abot(q)
else if (k > nzv) then
  kb=k-nzv
  if (kb < 1 .or. kb > ntop) then
    call ls_record_interp_oob()
    value=BOGUS
    return
  end if
  q=i+(j-1)*ldim+(kb-1)*ldim*nyv
  value=atop(q)
else
  q=i+(j-1)*ldim+(k-albzv)*ldim*nyv
  value=a(q)
end if
end function ls_halo_value

!*******************************************************************************
real(rprec) function ls_interp_field_halo(a,abot,atop,nbot,ntop,albzv,     &
    x,y,z,w_node,ldim,nxv,nyv,nzv,dxv,dyv,dzv,lxv,lyv) result(value)
!$acc routine seq
!*******************************************************************************
! Trilinear interpolation using the local field plus compact lower/upper-rank
! overlap arrays. z is expressed in the current rank's local coordinates.
real(rprec), intent(in) :: a(*), abot(*), atop(*)
integer, intent(in) :: nbot,ntop,albzv,ldim,nxv,nyv,nzv
real(rprec), intent(in) :: x,y,z,dxv,dyv,dzv,lxv,lyv
logical, intent(in) :: w_node
integer :: i0,i1,j0,j1,k0,k1
real(rprec) :: xm,ym,ax,ay,az,s,w00,w10,w01,w11
real(rprec) :: f000,f100,f010,f110,f001,f101,f011,f111

xm=modulo(x,lxv)
ym=modulo(y,lyv)
i0=min(nxv,max(1,int(floor(xm/dxv))+1))
j0=min(nyv,max(1,int(floor(ym/dyv))+1))
i1=modulo(i0,nxv)+1
j1=modulo(j0,nyv)+1
ax=modulo(xm,dxv)/dxv
ay=modulo(ym,dyv)/dyv
if (w_node) then
  s=1._rprec
else
  s=0.5_rprec
end if
k0=int(floor(z/dzv+s))
k1=k0+1
az=z/dzv-(floor(z/dzv+s)-s)

f000=ls_halo_value(a,abot,atop,nbot,ntop,albzv,k0,i0,j0,ldim,nyv,nzv)
f100=ls_halo_value(a,abot,atop,nbot,ntop,albzv,k0,i1,j0,ldim,nyv,nzv)
f010=ls_halo_value(a,abot,atop,nbot,ntop,albzv,k0,i0,j1,ldim,nyv,nzv)
f110=ls_halo_value(a,abot,atop,nbot,ntop,albzv,k0,i1,j1,ldim,nyv,nzv)
f001=ls_halo_value(a,abot,atop,nbot,ntop,albzv,k1,i0,j0,ldim,nyv,nzv)
f101=ls_halo_value(a,abot,atop,nbot,ntop,albzv,k1,i1,j0,ldim,nyv,nzv)
f011=ls_halo_value(a,abot,atop,nbot,ntop,albzv,k1,i0,j1,ldim,nyv,nzv)
f111=ls_halo_value(a,abot,atop,nbot,ntop,albzv,k1,i1,j1,ldim,nyv,nzv)
if (max(abs(f000),abs(f100),abs(f010),abs(f110),abs(f001),abs(f101),     &
        abs(f011),abs(f111)) >= 0.5_rprec*abs(BOGUS)) then
  value=BOGUS
  return
end if
w00=(1._rprec-ax)*(1._rprec-ay)
w10=ax*(1._rprec-ay)
w01=(1._rprec-ax)*ay
w11=ax*ay
value=(1._rprec-az)*(w00*f000+w10*f100+w01*f010+w11*f110) +          &
      az*(w00*f001+w10*f101+w01*f011+w11*f111)
end function ls_interp_field_halo

!*******************************************************************************
subroutine level_set_bc_gpu_core(phi_cutoff, phi_zero, normal)
!*******************************************************************************
real(rprec), intent(in) :: phi_cutoff, phi_zero
real(rprec), intent(in) :: normal(3,ld,ny,lbz:nz)

if (.not. use_log_profile) then
  if (use_extrap_tau_log) then
    call extrap_tau_log_gpu(phi_cutoff, normal)
  else
    call interp_tau_gpu(phi_cutoff, phi_zero, normal)
    if (use_extrap_tau_simple) then
      call extrap_tau_simple_gpu(phi_cutoff, phi_zero, normal)
    else
      call extrap_tau_legacy_gpu(phi_cutoff, phi_zero, normal)
    end if
  end if
end if
if (use_smooth_tau) then
  !$acc wait(1)
  call smooth_tau_gpu(-phi_cutoff)
end if
!$acc wait(1)
end subroutine level_set_bc_gpu_core

!$acc routine seq
!*******************************************************************************
subroutine modify_dutdn_gpu(i,j,k,tau,phix,xh,yh,zh,w_node,phi_cutoff,roughness,&
    kappa)
!*******************************************************************************
! Rotate the velocity-gradient tensor into the local wall frame, impose the
! log-law tangential/normal derivative, and rotate it back. This is the device
! equivalent of modify_dutdn; it is called only for the narrow interface band.
integer, intent(in) :: i,j,k
real(rprec), intent(in) :: tau,phix,xh(3),yh(3),zh(3),phi_cutoff,roughness,kappa
logical, intent(in) :: w_node
real(rprec) :: a(3,3),g(3,3),gp(3,3),tmp(3,3)
real(rprec) :: grad,phi_min
integer :: r,c,q

a(1,:)=xh
a(2,:)=yh
a(3,:)=zh
if (.not. w_node) then
  g(:,1)=(/dudx(i,j,k),dvdx(i,j,k),                                  &
           0.5_rprec*(dwdx(i,j,k)+dwdx(i,j,k+1))/)
  g(:,2)=(/dudy(i,j,k),dvdy(i,j,k),                                  &
           0.5_rprec*(dwdy(i,j,k)+dwdy(i,j,k+1))/)
  g(:,3)=(/0.5_rprec*(dudz(i,j,k)+dudz(i,j,k+1)),                     &
           0.5_rprec*(dvdz(i,j,k)+dvdz(i,j,k+1)),dwdz(i,j,k)/)
else
  g(:,1)=(/0.5_rprec*(dudx(i,j,k)+dudx(i,j,k-1)),                     &
           0.5_rprec*(dvdx(i,j,k)+dvdx(i,j,k-1)),dwdx(i,j,k)/)
  g(:,2)=(/0.5_rprec*(dudy(i,j,k)+dudy(i,j,k-1)),                     &
           0.5_rprec*(dvdy(i,j,k)+dvdy(i,j,k-1)),dwdy(i,j,k)/)
  g(:,3)=(/dudz(i,j,k),dvdz(i,j,k),                                  &
           0.5_rprec*(dwdz(i,j,k)+dwdz(i,j,k-1))/)
end if

do c=1,3
  do r=1,3
    tmp(r,c)=0._rprec
    do q=1,3
      tmp(r,c)=tmp(r,c)+g(r,q)*a(c,q)
    end do
  end do
end do
do c=1,3
  do r=1,3
    gp(r,c)=0._rprec
    do q=1,3
      gp(r,c)=gp(r,c)+a(r,q)*tmp(q,c)
    end do
  end do
end do

phi_min=0.1_rprec*phi_cutoff
grad=sqrt(abs(tau))/(kappa*(max(phi_min,phix)+roughness))
gp(1,3)=grad

do c=1,3
  do r=1,3
    tmp(r,c)=0._rprec
    do q=1,3
      tmp(r,c)=tmp(r,c)+gp(r,q)*a(q,c)
    end do
  end do
end do
do c=1,3
  do r=1,3
    g(r,c)=0._rprec
    do q=1,3
      g(r,c)=g(r,c)+a(q,r)*tmp(q,c)
    end do
  end do
end do

if (.not. w_node) then
  dudx(i,j,k)=g(1,1)
  dvdx(i,j,k)=g(2,1)
  dudy(i,j,k)=g(1,2)
  dvdy(i,j,k)=g(2,2)
  dwdz(i,j,k)=g(3,3)
else
  dudz(i,j,k)=g(1,3)
  dvdz(i,j,k)=g(2,3)
  dwdx(i,j,k)=g(3,1)
  dwdy(i,j,k)=g(3,2)
end if
end subroutine modify_dutdn_gpu

!*******************************************************************************
subroutine interp_tau_gpu(phi_cutoff, phi_zero, normal)
!*******************************************************************************
use param, only : vonK
real(rprec), intent(in) :: phi_cutoff, phi_zero
real(rprec), intent(in) :: normal(3,ld,ny,lbz:nz)
integer :: i, j, k, kmin
real(rprec), parameter :: eps = 100._rprec * epsilon(0._rprec)
real(rprec) :: phix, nxv, nyv, nzv, nmag, xv, yv, zv
  real(rprec) :: ux, uy, uz, dotn, txv, tyv, tzv, tmag, tau, wall_distance
  real(rprec) :: xh, yh, zh
  real(rprec) :: x_hat(3), y_hat(3), z_hat(3)

! uv-node stresses
!$acc parallel loop collapse(3) default(present) async(1)                    &
!$acc private(phix,nxv,nyv,nzv,nmag,xv,yv,zv,ux,uy,uz,dotn,txv,tyv,tzv,    &
!$acc         tmag,tau,wall_distance,xh,yh,zh,x_hat,y_hat,z_hat)
do k = 1, nz - 1
  do j = 1, ny
    do i = 1, nx
      phix = phi(i,j,k)
      if (phix >= phi_zero .and. phix <= phi_cutoff) then
        nxv = normal(1,i,j,k)
        nyv = normal(2,i,j,k)
        nzv = normal(3,i,j,k)
        nmag = sqrt(nxv*nxv + nyv*nyv + nzv*nzv)
        if (nmag > eps) then
          nxv = nxv/nmag
          nyv = nyv/nmag
          nzv = nzv/nmag
          if (physBC) then
            xv = real(i-1,rprec)*dx + nxv*(phi_cutoff-phix)
            yv = real(j-1,rprec)*dy + nyv*(phi_cutoff-phix)
            zv = (real(k,rprec)-0.5_rprec)*dz + nzv*(phi_cutoff-phix)
#ifdef PPMPI
            if (nproc > 1) then
              ux=ls_interp_field_halo(u,ubot,utop,nvelbot,nveltop,lbz,      &
                                      xv,yv,zv,.false.,ld,nx,ny,nz,         &
                                      dx,dy,dz,L_x,L_y)
              uy=ls_interp_field_halo(v,vbot,vtop,nvelbot,nveltop,lbz,      &
                                      xv,yv,zv,.false.,ld,nx,ny,nz,         &
                                      dx,dy,dz,L_x,L_y)
              uz=ls_interp_field_halo(w,wbot,wtop,nvelbot,nveltop,lbz,      &
                                      xv,yv,zv,.true.,ld,nx,ny,nz,          &
                                      dx,dy,dz,L_x,L_y)
            else
#endif
              ux=ls_interp_field(u,xv,yv,zv,.false.,LS_GRID_ARGS)
              uy=ls_interp_field(v,xv,yv,zv,.false.,LS_GRID_ARGS)
              uz=ls_interp_field(w,xv,yv,zv,.true.,LS_GRID_ARGS)
#ifdef PPMPI
            end if
#endif
          else
            ux = u(i,j,k)
            uy = v(i,j,k)
            uz = 0.5_rprec*(w(i,j,k)+w(i,j,k+1))
          end if
          dotn = ux*nxv + uy*nyv + uz*nzv
          txv = ux-dotn*nxv
          tyv = uy-dotn*nyv
          tzv = uz-dotn*nzv
          tmag = sqrt(txv*txv + tyv*tyv + tzv*tzv)
          if (tmag <= eps) then
            txx(i,j,k)=0._rprec
            txy(i,j,k)=0._rprec
            tyy(i,j,k)=0._rprec
            tzz(i,j,k)=0._rprec
          else
            xh=txv/tmag
            yh=tyv/tmag
            zh=tzv/tmag
            if (physBC) then
              tau = -(vonK*tmag/log(1._rprec+phi_cutoff/zo_level_set))**2
            else
              wall_distance=max(phix,zo_level_set)
              tau = -(vonK*tmag/log(1._rprec+wall_distance/zo_level_set))**2
            end if
            txx(i,j,k)=2._rprec*xh*nxv*tau
            txy(i,j,k)=(xh*nyv+nxv*yh)*tau
            tyy(i,j,k)=2._rprec*yh*nyv*tau
            tzz(i,j,k)=2._rprec*zh*nzv*tau
            if (use_modify_dutdn) then
              x_hat=(/xh,yh,zh/)
              z_hat=(/nxv,nyv,nzv/)
              y_hat=(/nyv*zh-nzv*yh,nzv*xh-nxv*zh,nxv*yh-nyv*xh/)
              call modify_dutdn_gpu(i,j,k,tau,phix,x_hat,y_hat,z_hat,       &
                                    .false.,phi_cutoff,zo_level_set,vonK)
            end if
          end if
        end if
      end if
    end do
  end do
end do

if (coord == 0) then
  kmin = 2
else
  kmin = 1
end if
! w-node stresses
!$acc parallel loop collapse(3) default(present) async(1)                    &
!$acc private(phix,nxv,nyv,nzv,nmag,xv,yv,zv,ux,uy,uz,dotn,txv,tyv,tzv,    &
!$acc         tmag,tau,wall_distance,xh,yh,zh,x_hat,y_hat,z_hat)
do k = kmin, nz - 1
  do j = 1, ny
    do i = 1, nx
      phix = 0.5_rprec*(phi(i,j,k)+phi(i,j,k-1))
      if (phix >= phi_zero .and. phix <= phi_cutoff) then
        nxv=0.5_rprec*(normal(1,i,j,k)+normal(1,i,j,k-1))
        nyv=0.5_rprec*(normal(2,i,j,k)+normal(2,i,j,k-1))
        nzv=0.5_rprec*(normal(3,i,j,k)+normal(3,i,j,k-1))
        nmag=sqrt(nxv*nxv+nyv*nyv+nzv*nzv)
        if (nmag > eps) then
          nxv=nxv/nmag
          nyv=nyv/nmag
          nzv=nzv/nmag
          if (physBC) then
            xv=real(i-1,rprec)*dx+nxv*(phi_cutoff-phix)
            yv=real(j-1,rprec)*dy+nyv*(phi_cutoff-phix)
            zv=real(k-1,rprec)*dz+nzv*(phi_cutoff-phix)
#ifdef PPMPI
            if (nproc > 1) then
              ux=ls_interp_field_halo(u,ubot,utop,nvelbot,nveltop,lbz,      &
                                      xv,yv,zv,.false.,ld,nx,ny,nz,         &
                                      dx,dy,dz,L_x,L_y)
              uy=ls_interp_field_halo(v,vbot,vtop,nvelbot,nveltop,lbz,      &
                                      xv,yv,zv,.false.,ld,nx,ny,nz,         &
                                      dx,dy,dz,L_x,L_y)
              uz=ls_interp_field_halo(w,wbot,wtop,nvelbot,nveltop,lbz,      &
                                      xv,yv,zv,.true.,ld,nx,ny,nz,          &
                                      dx,dy,dz,L_x,L_y)
            else
#endif
              ux=ls_interp_field(u,xv,yv,zv,.false.,LS_GRID_ARGS)
              uy=ls_interp_field(v,xv,yv,zv,.false.,LS_GRID_ARGS)
              uz=ls_interp_field(w,xv,yv,zv,.true.,LS_GRID_ARGS)
#ifdef PPMPI
            end if
#endif
          else
            ux=0.5_rprec*(u(i,j,k)+u(i,j,k-1))
            uy=0.5_rprec*(v(i,j,k)+v(i,j,k-1))
            uz=w(i,j,k)
          end if
          dotn=ux*nxv+uy*nyv+uz*nzv
          txv=ux-dotn*nxv
          tyv=uy-dotn*nyv
          tzv=uz-dotn*nzv
          tmag=sqrt(txv*txv+tyv*tyv+tzv*tzv)
          if (tmag <= eps) then
            txz(i,j,k)=0._rprec
            tyz(i,j,k)=0._rprec
          else
            xh=txv/tmag
            yh=tyv/tmag
            zh=tzv/tmag
            if (physBC) then
              tau=-(vonK*tmag/log(1._rprec+phi_cutoff/zo_level_set))**2
            else
              wall_distance=max(phix,zo_level_set)
              tau=-(vonK*tmag/log(1._rprec+wall_distance/zo_level_set))**2
            end if
            txz(i,j,k)=(xh*nzv+nxv*zh)*tau
            tyz(i,j,k)=(yh*nzv+nyv*zh)*tau
            if (use_modify_dutdn) then
              x_hat=(/xh,yh,zh/)
              z_hat=(/nxv,nyv,nzv/)
              y_hat=(/nyv*zh-nzv*yh,nzv*xh-nxv*zh,nxv*yh-nyv*xh/)
              call modify_dutdn_gpu(i,j,k,tau,phix,x_hat,y_hat,z_hat,       &
                                    .true.,phi_cutoff,zo_level_set,vonK)
            end if
          end if
        end if
      end if
    end do
  end do
end do
end subroutine interp_tau_gpu

!*******************************************************************************
subroutine level_set_desired_velocity_gpu(phi_zero,phi_cutoff,normal,          &
    do_enforce_un,do_log_profile,udes,vdes,wdes)
!*******************************************************************************
! Construct the optional immersed-boundary target velocity entirely on device.
! If both historical switches are enabled, the log-profile result wins, just as
! in the CPU sequence where enforce_log_profile is called second and reinitializes
! all three target arrays.
use test_filtermodule, only : filter_size
real(rprec), intent(in) :: phi_zero,phi_cutoff
real(rprec), intent(in) :: normal(3,ld,ny,lbz:nz)
logical, intent(in) :: do_enforce_un,do_log_profile
real(rprec), intent(inout) :: udes(ld,ny,lbz:nz),vdes(ld,ny,lbz:nz),           &
                              wdes(ld,ny,lbz:nz)
integer :: i,j,k
real(rprec) :: dphi,phi1,phi2,x1,y1,z1,x2,y2,z2
real(rprec) :: nxv,nyv,nzv,nmag,ux,uy,uz,dotn,txv,tyv,tzv
real(rprec) :: normal_ratio,tangent_ratio

!$acc parallel loop collapse(3) default(present) async(1)
do k=lbz,nz
  do j=1,ny
    do i=1,ld
      udes(i,j,k)=huge(1._rprec)
      vdes(i,j,k)=huge(1._rprec)
      wdes(i,j,k)=huge(1._rprec)
    end do
  end do
end do
if (.not. do_enforce_un .and. .not. do_log_profile) return

dphi=filter_size*sqrt(dx*dx+dy*dy+dz*dz)
! uv-node target components
!$acc parallel loop collapse(3) default(present) async(1)                    &
!$acc private(phi1,phi2,x1,y1,z1,x2,y2,z2,nxv,nyv,nzv,ux,uy,uz,dotn,       &
!$acc         txv,tyv,tzv,normal_ratio,tangent_ratio)
do k=1,nz-1
  do j=1,ny
    do i=1,nx
      phi1=phi(i,j,k)
      if (phi_zero < phi1 .and. phi1 < phi_cutoff) then
        x1=real(i-1,rprec)*dx
        y1=real(j-1,rprec)*dy
        z1=(real(k,rprec)-0.5_rprec)*dz
        nxv=normal(1,i,j,k)
        nyv=normal(2,i,j,k)
        nzv=normal(3,i,j,k)
        if (do_log_profile) then
          x2=x1+dphi*nxv
          y2=y1+dphi*nyv
          z2=z1+dphi*nzv
          phi2=phi1+dphi
#ifdef PPMPI
          if (nproc > 1) then
            ux=ls_interp_field_halo(u,ubot,utop,nvelbot,nveltop,lbz,         &
                                    x2,y2,z2,.false.,LS_HALO_GRID_ARGS)
            uy=ls_interp_field_halo(v,vbot,vtop,nvelbot,nveltop,lbz,         &
                                    x2,y2,z2,.false.,LS_HALO_GRID_ARGS)
            uz=ls_interp_field_halo(w,wbot,wtop,nvelbot,nveltop,lbz,         &
                                    x2,y2,z2,.true.,LS_HALO_GRID_ARGS)
          else
#endif
            ux=ls_interp_field(u,x2,y2,z2,.false.,LS_GRID_ARGS)
            uy=ls_interp_field(v,x2,y2,z2,.false.,LS_GRID_ARGS)
            uz=ls_interp_field(w,x2,y2,z2,.true.,LS_GRID_ARGS)
#ifdef PPMPI
          end if
#endif
          dotn=ux*nxv+uy*nyv+uz*nzv
          txv=ux-dotn*nxv
          tyv=uy-dotn*nyv
          tzv=uz-dotn*nzv
          normal_ratio=(phi1/phi2)**2
          tangent_ratio=log(1._rprec+phi1/zo_level_set) /                    &
                        log(1._rprec+phi2/zo_level_set)
          udes(i,j,k)=tangent_ratio*txv+normal_ratio*dotn*nxv
          vdes(i,j,k)=tangent_ratio*tyv+normal_ratio*dotn*nyv
        else
          x2=x1+dphi*nxv
          y2=y1+dphi*nyv
          z2=z1+dphi*nzv
          phi2=phi1+dphi
#ifdef PPMPI
          if (nproc > 1) then
            txv=ls_interp_field_halo(u,ubot,utop,nvelbot,nveltop,lbz,        &
                                     x1,y1,z1,.false.,LS_HALO_GRID_ARGS)
            tyv=ls_interp_field_halo(v,vbot,vtop,nvelbot,nveltop,lbz,        &
                                     x1,y1,z1,.false.,LS_HALO_GRID_ARGS)
            tzv=ls_interp_field_halo(w,wbot,wtop,nvelbot,nveltop,lbz,        &
                                     x1,y1,z1,.true.,LS_HALO_GRID_ARGS)
            ux=ls_interp_field_halo(u,ubot,utop,nvelbot,nveltop,lbz,         &
                                    x2,y2,z2,.false.,LS_HALO_GRID_ARGS)
            uy=ls_interp_field_halo(v,vbot,vtop,nvelbot,nveltop,lbz,         &
                                    x2,y2,z2,.false.,LS_HALO_GRID_ARGS)
            uz=ls_interp_field_halo(w,wbot,wtop,nvelbot,nveltop,lbz,         &
                                    x2,y2,z2,.true.,LS_HALO_GRID_ARGS)
          else
#endif
            txv=ls_interp_field(u,x1,y1,z1,.false.,LS_GRID_ARGS)
            tyv=ls_interp_field(v,x1,y1,z1,.false.,LS_GRID_ARGS)
            tzv=ls_interp_field(w,x1,y1,z1,.true.,LS_GRID_ARGS)
            ux=ls_interp_field(u,x2,y2,z2,.false.,LS_GRID_ARGS)
            uy=ls_interp_field(v,x2,y2,z2,.false.,LS_GRID_ARGS)
            uz=ls_interp_field(w,x2,y2,z2,.true.,LS_GRID_ARGS)
#ifdef PPMPI
          end if
#endif
          dotn=ux*nxv+uy*nyv+uz*nzv
          tangent_ratio=txv*nxv+tyv*nyv+tzv*nzv
          normal_ratio=(phi1/phi2)**2
          udes(i,j,k)=txv+(normal_ratio*dotn-tangent_ratio)*nxv
          vdes(i,j,k)=tyv+(normal_ratio*dotn-tangent_ratio)*nyv
        end if
      end if
    end do
  end do
end do

! w-node target component
!$acc parallel loop collapse(3) default(present) async(1)                    &
!$acc private(phi1,phi2,x1,y1,z1,x2,y2,z2,nxv,nyv,nzv,nmag,ux,uy,uz,       &
!$acc         dotn,txv,tyv,tzv,normal_ratio,tangent_ratio)
do k=2,nz-1
  do j=1,ny
    do i=1,nx
      phi1=0.5_rprec*(phi(i,j,k)+phi(i,j,k-1))
      if (phi_zero < phi1 .and. phi1 < phi_cutoff) then
        x1=real(i-1,rprec)*dx
        y1=real(j-1,rprec)*dy
        z1=real(k-1,rprec)*dz
        nxv=0.5_rprec*(normal(1,i,j,k)+normal(1,i,j,k-1))
        nyv=0.5_rprec*(normal(2,i,j,k)+normal(2,i,j,k-1))
        nzv=0.5_rprec*(normal(3,i,j,k)+normal(3,i,j,k-1))
        nmag=sqrt(nxv*nxv+nyv*nyv+nzv*nzv)
        if (nmag > 0._rprec) then
          nxv=nxv/nmag
          nyv=nyv/nmag
          nzv=nzv/nmag
        end if
        if (do_log_profile) then
          x2=x1+dphi*nxv
          y2=y1+dphi*nyv
          z2=z1+dphi*nzv
          phi2=phi1+dphi
#ifdef PPMPI
          if (nproc > 1) then
            ux=ls_interp_field_halo(u,ubot,utop,nvelbot,nveltop,lbz,         &
                                    x2,y2,z2,.false.,LS_HALO_GRID_ARGS)
            uy=ls_interp_field_halo(v,vbot,vtop,nvelbot,nveltop,lbz,         &
                                    x2,y2,z2,.false.,LS_HALO_GRID_ARGS)
            uz=ls_interp_field_halo(w,wbot,wtop,nvelbot,nveltop,lbz,         &
                                    x2,y2,z2,.true.,LS_HALO_GRID_ARGS)
          else
#endif
            ux=ls_interp_field(u,x2,y2,z2,.false.,LS_GRID_ARGS)
            uy=ls_interp_field(v,x2,y2,z2,.false.,LS_GRID_ARGS)
            uz=ls_interp_field(w,x2,y2,z2,.true.,LS_GRID_ARGS)
#ifdef PPMPI
          end if
#endif
          dotn=ux*nxv+uy*nyv+uz*nzv
          txv=ux-dotn*nxv
          tyv=uy-dotn*nyv
          tzv=uz-dotn*nzv
          normal_ratio=(phi1/phi2)**2
          tangent_ratio=log(1._rprec+phi1/zo_level_set) /                    &
                        log(1._rprec+phi2/zo_level_set)
          wdes(i,j,k)=tangent_ratio*tzv+normal_ratio*dotn*nzv
        else
          x2=x1+dphi*nxv
          y2=y1+dphi*nyv
          z2=z1+dphi*nzv
          phi2=phi1+dphi
#ifdef PPMPI
          if (nproc > 1) then
            txv=ls_interp_field_halo(u,ubot,utop,nvelbot,nveltop,lbz,        &
                                     x1,y1,z1,.false.,LS_HALO_GRID_ARGS)
            tyv=ls_interp_field_halo(v,vbot,vtop,nvelbot,nveltop,lbz,        &
                                     x1,y1,z1,.false.,LS_HALO_GRID_ARGS)
            tzv=ls_interp_field_halo(w,wbot,wtop,nvelbot,nveltop,lbz,        &
                                     x1,y1,z1,.true.,LS_HALO_GRID_ARGS)
            ux=ls_interp_field_halo(u,ubot,utop,nvelbot,nveltop,lbz,         &
                                    x2,y2,z2,.false.,LS_HALO_GRID_ARGS)
            uy=ls_interp_field_halo(v,vbot,vtop,nvelbot,nveltop,lbz,         &
                                    x2,y2,z2,.false.,LS_HALO_GRID_ARGS)
            uz=ls_interp_field_halo(w,wbot,wtop,nvelbot,nveltop,lbz,         &
                                    x2,y2,z2,.true.,LS_HALO_GRID_ARGS)
          else
#endif
            txv=ls_interp_field(u,x1,y1,z1,.false.,LS_GRID_ARGS)
            tyv=ls_interp_field(v,x1,y1,z1,.false.,LS_GRID_ARGS)
            tzv=ls_interp_field(w,x1,y1,z1,.true.,LS_GRID_ARGS)
            ux=ls_interp_field(u,x2,y2,z2,.false.,LS_GRID_ARGS)
            uy=ls_interp_field(v,x2,y2,z2,.false.,LS_GRID_ARGS)
            uz=ls_interp_field(w,x2,y2,z2,.true.,LS_GRID_ARGS)
#ifdef PPMPI
          end if
#endif
          dotn=ux*nxv+uy*nyv+uz*nzv
          tangent_ratio=txv*nxv+tyv*nyv+tzv*nzv
          normal_ratio=(phi1/phi2)**2
          wdes(i,j,k)=tzv+(normal_ratio*dotn-tangent_ratio)*nzv
        end if
      end if
    end do
  end do
end do
end subroutine level_set_desired_velocity_gpu

!*******************************************************************************
subroutine extrap_tau_log_gpu(phi_cutoff,normal)
!*******************************************************************************
! Log-law image-point stress treatment. The CPU implementation does not support
! cross-rank interpolation for this option, so the caller restricts it to one
! rank; all interpolation and extrapolation still remain device resident.
use param, only : vonK
real(rprec), intent(in) :: phi_cutoff
real(rprec), intent(in) :: normal(3,ld,ny,lbz:nz)
integer :: i,j,k
real(rprec), parameter :: eps=100._rprec*epsilon(1._rprec)
real(rprec) :: phi_a,phix,phi1,phi2,dphi,x,y,z,x1,y1,z1,x2,y2,z2
real(rprec) :: nxv,nyv,nzv,nmag,ux,uy,uz,vn,txv,tyv,tzv,tmag
real(rprec) :: xh,yh,zh,tau_w,wgt,wgt_im
real(rprec) :: a1,b1,c1,d1,a2,b2,c2,d2,aim,bim,cim,dim

phi_a=-phi_cutoff
dphi=phi_cutoff

! uv-node components: txx, txy, tyy, tzz
!$acc parallel loop collapse(3) default(present) async(1)                    &
!$acc private(phix,phi1,phi2,x,y,z,x1,y1,z1,x2,y2,z2,nxv,nyv,nzv,         &
!$acc         ux,uy,uz,vn,txv,tyv,tzv,tmag,xh,yh,zh,tau_w,wgt,wgt_im,     &
!$acc         a1,b1,c1,d1,a2,b2,c2,d2,aim,bim,cim,dim)
do k=1,nz-1
  do j=1,ny
    do i=1,nx
      phix=phi(i,j,k)
      if (phi_a < phix .and. phix < 0._rprec) then
        x=real(i-1,rprec)*dx
        y=real(j-1,rprec)*dy
        z=(real(k,rprec)-0.5_rprec)*dz
        nxv=normal(1,i,j,k)
        nyv=normal(2,i,j,k)
        nzv=normal(3,i,j,k)
        x1=x+(dphi-phi_a)*nxv
        y1=y+(dphi-phi_a)*nyv
        z1=z+(dphi-phi_a)*nzv
        ux=ls_interp_field(u,x1,y1,z1,.false.,LS_GRID_ARGS)
        uy=ls_interp_field(v,x1,y1,z1,.false.,LS_GRID_ARGS)
        uz=ls_interp_field(w,x1,y1,z1,.true., LS_GRID_ARGS)
        vn=ux*nxv+uy*nyv+uz*nzv
        txv=ux-vn*nxv
        tyv=uy-vn*nyv
        tzv=uz-vn*nzv
        tmag=sqrt(txv*txv+tyv*tyv+tzv*tzv)
        if (tmag < eps) then
          txx(i,j,k)=0._rprec
          txy(i,j,k)=0._rprec
          tyy(i,j,k)=0._rprec
          tzz(i,j,k)=0._rprec
        else
          xh=txv/tmag
          yh=tyv/tmag
          zh=tzv/tmag
          phi1=ls_interp_field(phi,x1,y1,z1,.false.,LS_GRID_ARGS)
          tau_w=-(tmag*vonK/log(phi1/zo_level_set))**2
          a1=ls_interp_field(txx,x1,y1,z1,.false.,LS_GRID_ARGS)
          b1=ls_interp_field(txy,x1,y1,z1,.false.,LS_GRID_ARGS)
          c1=ls_interp_field(tyy,x1,y1,z1,.false.,LS_GRID_ARGS)
          d1=ls_interp_field(tzz,x1,y1,z1,.false.,LS_GRID_ARGS)
          wgt=abs(phi1)/(abs(phix)+abs(phi1))
          if (wgt >= 0.5_rprec) then
            txx(i,j,k)=(2._rprec*xh*nxv*tau_w-(1._rprec-wgt)*a1)/wgt
            txy(i,j,k)=((xh*nyv+nxv*yh)*tau_w-(1._rprec-wgt)*b1)/wgt
            tyy(i,j,k)=(2._rprec*yh*nyv*tau_w-(1._rprec-wgt)*c1)/wgt
            tzz(i,j,k)=(2._rprec*zh*nzv*tau_w-(1._rprec-wgt)*d1)/wgt
          else
            x2=x1+dphi*nxv
            y2=y1+dphi*nyv
            z2=z1+dphi*nzv
            phi2=ls_interp_field(phi,x2,y2,z2,.false.,LS_GRID_ARGS)
            a2=ls_interp_field(txx,x2,y2,z2,.false.,LS_GRID_ARGS)
            b2=ls_interp_field(txy,x2,y2,z2,.false.,LS_GRID_ARGS)
            c2=ls_interp_field(tyy,x2,y2,z2,.false.,LS_GRID_ARGS)
            d2=ls_interp_field(tzz,x2,y2,z2,.false.,LS_GRID_ARGS)
            wgt_im=(abs(phi2)-abs(phix))/(abs(phi2)-abs(phi1))
            aim=wgt_im*a1+(1._rprec-wgt_im)*a2
            bim=wgt_im*b1+(1._rprec-wgt_im)*b2
            cim=wgt_im*c1+(1._rprec-wgt_im)*c2
            dim=wgt_im*d1+(1._rprec-wgt_im)*d2
            txx(i,j,k)=4._rprec*xh*nxv*tau_w-aim
            txy(i,j,k)=2._rprec*(xh*nyv+nxv*yh)*tau_w-bim
            tyy(i,j,k)=4._rprec*yh*nyv*tau_w-cim
            tzz(i,j,k)=4._rprec*zh*nzv*tau_w-dim
          end if
        end if
      end if
    end do
  end do
end do

! w-node components: txz, tyz
!$acc parallel loop collapse(3) default(present) async(1)                    &
!$acc private(phix,phi1,phi2,x,y,z,x1,y1,z1,x2,y2,z2,nxv,nyv,nzv,nmag,    &
!$acc         ux,uy,uz,vn,txv,tyv,tzv,tmag,xh,yh,zh,tau_w,wgt,wgt_im,     &
!$acc         a1,b1,a2,b2,aim,bim)
do k=2,nz
  do j=1,ny
    do i=1,nx
      phix=0.5_rprec*(phi(i,j,k)+phi(i,j,k-1))
      if (phi_a < phix .and. phix < 0._rprec) then
        x=real(i-1,rprec)*dx
        y=real(j-1,rprec)*dy
        z=real(k-1,rprec)*dz
        nxv=0.5_rprec*(normal(1,i,j,k)+normal(1,i,j,k-1))
        nyv=0.5_rprec*(normal(2,i,j,k)+normal(2,i,j,k-1))
        nzv=0.5_rprec*(normal(3,i,j,k)+normal(3,i,j,k-1))
        nmag=sqrt(nxv*nxv+nyv*nyv+nzv*nzv)
        if (nmag > 0._rprec) then
          nxv=nxv/nmag
          nyv=nyv/nmag
          nzv=nzv/nmag
        end if
        x1=x+(dphi-phi_a)*nxv
        y1=y+(dphi-phi_a)*nyv
        z1=z+(dphi-phi_a)*nzv
        ux=ls_interp_field(u,x1,y1,z1,.false.,LS_GRID_ARGS)
        uy=ls_interp_field(v,x1,y1,z1,.false.,LS_GRID_ARGS)
        uz=ls_interp_field(w,x1,y1,z1,.true., LS_GRID_ARGS)
        vn=ux*nxv+uy*nyv+uz*nzv
        txv=ux-vn*nxv
        tyv=uy-vn*nyv
        tzv=uz-vn*nzv
        tmag=sqrt(txv*txv+tyv*tyv+tzv*tzv)
        if (tmag < eps) then
          txz(i,j,k)=0._rprec
          tyz(i,j,k)=0._rprec
        else
          xh=txv/tmag
          yh=tyv/tmag
          zh=tzv/tmag
          phi1=ls_interp_field(phi,x1,y1,z1,.false.,LS_GRID_ARGS)
          tau_w=-(tmag*vonK/log(phi1/zo_level_set))**2
          a1=ls_interp_field(txz,x1,y1,z1,.true.,LS_GRID_ARGS)
          b1=ls_interp_field(tyz,x1,y1,z1,.true.,LS_GRID_ARGS)
          wgt=abs(phi1)/(abs(phix)+abs(phi1))
          if (wgt >= 0.5_rprec) then
            txz(i,j,k)=((xh*nzv+nxv*zh)*tau_w-(1._rprec-wgt)*a1)/wgt
            tyz(i,j,k)=((yh*nzv+nyv*zh)*tau_w-(1._rprec-wgt)*b1)/wgt
          else
            x2=x1+dphi*nxv
            y2=y1+dphi*nyv
            z2=z1+dphi*nzv
            phi2=ls_interp_field(phi,x2,y2,z2,.false.,LS_GRID_ARGS)
            a2=ls_interp_field(txz,x2,y2,z2,.true.,LS_GRID_ARGS)
            b2=ls_interp_field(tyz,x2,y2,z2,.true.,LS_GRID_ARGS)
            wgt_im=(abs(phi2)-abs(phix))/(abs(phi2)-abs(phi1))
            aim=wgt_im*a1+(1._rprec-wgt_im)*a2
            bim=wgt_im*b1+(1._rprec-wgt_im)*b2
            txz(i,j,k)=2._rprec*(xh*nzv+nxv*zh)*tau_w-aim
            tyz(i,j,k)=2._rprec*(yh*nzv+nyv*zh)*tau_w-bim
          end if
        end if
      end if
    end do
  end do
end do
end subroutine extrap_tau_log_gpu

!$acc routine seq
!*******************************************************************************
real(rprec) function ls_fit3_value(a,pi,pj,pk,li,lj,lk,nlist,ldim,nyv,lbzv, &
    dxv,dyv,dzv) result(value)
!*******************************************************************************
! Evaluate the legacy local plane fit at the target point. The fit uses the
! first three fluid corners lying in one coordinate plane, matching fit3().
real(rprec), intent(in) :: a(*)
integer, intent(in) :: pi,pj,pk,nlist,li(7),lj(7),lk(7)
integer, intent(in) :: ldim,nyv,lbzv
real(rprec), intent(in) :: dxv,dyv,dzv
integer :: m,d,count,dir,q,sel(3),i1,i2,i3
real(rprec) :: x1,x2,x3,y1,y2,y3,t1,t2,t3,det,det0,total

if (nlist == 1) then
  q=li(1)+(lj(1)-1)*ldim+(lk(1)-lbzv)*ldim*nyv
  value=a(q)
  return
else if (nlist == 2) then
  q=li(1)+(lj(1)-1)*ldim+(lk(1)-lbzv)*ldim*nyv
  value=a(q)
  q=li(2)+(lj(2)-1)*ldim+(lk(2)-lbzv)*ldim*nyv
  value=0.5_rprec*(value+a(q))
  return
end if

dir=0
do d=1,3
  count=0
  do m=1,nlist
    if ((d == 1 .and. li(m) == pi) .or.                               &
        (d == 2 .and. lj(m) == pj) .or.                               &
        (d == 3 .and. lk(m) == pk)) count=count+1
  end do
  if (count >= 3) then
    dir=d
    exit
  end if
end do

if (dir == 0) then
  total=0._rprec
  do m=1,nlist
    q=li(m)+(lj(m)-1)*ldim+(lk(m)-lbzv)*ldim*nyv
    total=total+a(q)
  end do
  value=total/real(nlist,rprec)
  return
end if

count=0
do m=1,nlist
  if ((dir == 1 .and. li(m) == pi) .or.                               &
      (dir == 2 .and. lj(m) == pj) .or.                               &
      (dir == 3 .and. lk(m) == pk)) then
    count=count+1
    sel(count)=m
    if (count == 3) exit
  end if
end do
i1=sel(1)
i2=sel(2)
i3=sel(3)
select case(dir)
case(1)
  x1=real(lj(i1)-pj,rprec)*dyv; y1=real(lk(i1)-pk,rprec)*dzv
  x2=real(lj(i2)-pj,rprec)*dyv; y2=real(lk(i2)-pk,rprec)*dzv
  x3=real(lj(i3)-pj,rprec)*dyv; y3=real(lk(i3)-pk,rprec)*dzv
case(2)
  x1=real(li(i1)-pi,rprec)*dxv; y1=real(lk(i1)-pk,rprec)*dzv
  x2=real(li(i2)-pi,rprec)*dxv; y2=real(lk(i2)-pk,rprec)*dzv
  x3=real(li(i3)-pi,rprec)*dxv; y3=real(lk(i3)-pk,rprec)*dzv
case default
  x1=real(li(i1)-pi,rprec)*dxv; y1=real(lj(i1)-pj,rprec)*dyv
  x2=real(li(i2)-pi,rprec)*dxv; y2=real(lj(i2)-pj,rprec)*dyv
  x3=real(li(i3)-pi,rprec)*dxv; y3=real(lj(i3)-pj,rprec)*dyv
end select
q=li(i1)+(lj(i1)-1)*ldim+(lk(i1)-lbzv)*ldim*nyv; t1=a(q)
q=li(i2)+(lj(i2)-1)*ldim+(lk(i2)-lbzv)*ldim*nyv; t2=a(q)
q=li(i3)+(lj(i3)-1)*ldim+(lk(i3)-lbzv)*ldim*nyv; t3=a(q)
det=x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)
if (abs(det) <= tiny(1._rprec)) then
  value=(t1+t2+t3)/3._rprec
else
  det0=t1*(x2*y3-y2*x3)-x1*(t2*y3-y2*t3)+y1*(t2*x3-x2*t3)
  value=det0/det
end if
end function ls_fit3_value

!*******************************************************************************
subroutine extrap_tau_legacy_gpu(phi_cutoff,phi_zero,normal)
!*******************************************************************************
! Device implementation of the original seven-corner extrapolator. Candidate
! values always come from fluid nodes, so target points are independent and no
! frozen full-domain stress copy is required.
real(rprec), intent(in) :: phi_cutoff,phi_zero
real(rprec), intent(in) :: normal(3,ld,ny,lbz:nz)
integer :: i,j,k,m,di,dj,dk,si,sj,sk,ii,jj,kk,nlist
integer :: li(7),lj(7),lk(7)
real(rprec) :: phiw,nxv,nyv,nzv

! uv-node stresses
!$acc parallel loop collapse(3) default(present) async(1)                    &
!$acc private(m,di,dj,dk,si,sj,sk,ii,jj,kk,nlist,li,lj,lk)
do k=2,nz-1
  do j=2,ny-1
    do i=2,nx-1
      if (phi(i,j,k) < phi_zero .and. phi(i,j,k) >= -phi_cutoff) then
        si=1; if (normal(1,i,j,k) < 0._rprec) si=-1
        sj=1; if (normal(2,i,j,k) < 0._rprec) sj=-1
        sk=1; if (normal(3,i,j,k) < 0._rprec) sk=-1
        nlist=0
        do m=1,7
          di=modulo(m,2)
          dj=modulo(m/2,2)
          dk=modulo(m/4,2)
          ii=i+di*si
          jj=j+dj*sj
          kk=k+dk*sk
          if (phi(ii,jj,kk) >= 0._rprec) then
            nlist=nlist+1
            li(nlist)=ii; lj(nlist)=jj; lk(nlist)=kk
          end if
        end do
        if (nlist > 0) then
          txx(i,j,k)=ls_fit3_value(txx,i,j,k,li,lj,lk,nlist,ld,ny,lbz,dx,dy,dz)
          txy(i,j,k)=ls_fit3_value(txy,i,j,k,li,lj,lk,nlist,ld,ny,lbz,dx,dy,dz)
          tyy(i,j,k)=ls_fit3_value(tyy,i,j,k,li,lj,lk,nlist,ld,ny,lbz,dx,dy,dz)
          tzz(i,j,k)=ls_fit3_value(tzz,i,j,k,li,lj,lk,nlist,ld,ny,lbz,dx,dy,dz)
        end if
      end if
    end do
  end do
end do

! w-node stresses
!$acc parallel loop collapse(3) default(present) async(1)                    &
!$acc private(m,di,dj,dk,si,sj,sk,ii,jj,kk,nlist,li,lj,lk,phiw,nxv,nyv,nzv)
do k=2,nz-1
  do j=2,ny-1
    do i=2,nx-1
      phiw=0.5_rprec*(phi(i,j,k)+phi(i,j,k-1))
      if (phiw < phi_zero .and. phiw >= -phi_cutoff) then
        nxv=0.5_rprec*(normal(1,i,j,k)+normal(1,i,j,k-1))
        nyv=0.5_rprec*(normal(2,i,j,k)+normal(2,i,j,k-1))
        nzv=0.5_rprec*(normal(3,i,j,k)+normal(3,i,j,k-1))
        si=1; if (nxv < 0._rprec) si=-1
        sj=1; if (nyv < 0._rprec) sj=-1
        sk=1; if (nzv < 0._rprec) sk=-1
        nlist=0
        do m=1,7
          di=modulo(m,2)
          dj=modulo(m/2,2)
          dk=modulo(m/4,2)
          ii=i+di*si
          jj=j+dj*sj
          kk=k+dk*sk
          phiw=0.5_rprec*(phi(ii,jj,kk)+phi(ii,jj,kk-1))
          if (phiw >= 0._rprec) then
            nlist=nlist+1
            li(nlist)=ii; lj(nlist)=jj; lk(nlist)=kk
          end if
        end do
        if (nlist > 0) then
          txz(i,j,k)=ls_fit3_value(txz,i,j,k,li,lj,lk,nlist,ld,ny,lbz,dx,dy,dz)
          tyz(i,j,k)=ls_fit3_value(tyz,i,j,k,li,lj,lk,nlist,ld,ny,lbz,dx,dy,dz)
        end if
      end if
    end do
  end do
end do
end subroutine extrap_tau_legacy_gpu

!*******************************************************************************
subroutine extrap_tau_simple_gpu(phi_cutoff, phi_zero, normal)
!*******************************************************************************
real(rprec), intent(in) :: phi_cutoff, phi_zero
real(rprec), intent(in) :: normal(3,ld,ny,lbz:nz)
call extrap_tau_field_gpu(txx,txxbot,txxtop,.false.,phi_cutoff,phi_zero,normal)
call extrap_tau_field_gpu(txy,txybot,txytop,.false.,phi_cutoff,phi_zero,normal)
call extrap_tau_field_gpu(tyy,tyybot,tyytop,.false.,phi_cutoff,phi_zero,normal)
call extrap_tau_field_gpu(tzz,tzzbot,tzztop,.false.,phi_cutoff,phi_zero,normal)
call extrap_tau_field_gpu(txz,txzbot,txztop,.true., phi_cutoff,phi_zero,normal)
call extrap_tau_field_gpu(tyz,tyzbot,tyztop,.true., phi_cutoff,phi_zero,normal)
end subroutine extrap_tau_simple_gpu

!*******************************************************************************
subroutine extrap_tau_field_gpu(a,abot,atop,w_node,phi_cutoff,phi_zero,normal)
!*******************************************************************************
use test_filtermodule, only : filter_size
real(rprec), intent(inout) :: a(ld,ny,lbz:nz)
real(rprec), intent(in) :: abot(ld,ny,ntaubot), atop(ld,ny,ntautop)
logical, intent(in) :: w_node
real(rprec), intent(in) :: phi_cutoff, phi_zero
real(rprec), intent(in) :: normal(3,ld,ny,lbz:nz)
integer :: i,j,k,kmin
real(rprec), parameter :: eps=100._rprec*epsilon(0._rprec)
real(rprec) :: phix,nxv,nyv,nzv,nmag,x,y,z,dstep
real(rprec) :: x1,y1,z1,x2,y2,z2,phi1,a1,a2

! Freeze the complete source field before any immersed target is overwritten.
!$acc parallel loop collapse(3) present(a,tau_source) async(1)
do k=lbz,nz
  do j=1,ny
    do i=1,ld
      tau_source(i,j,k)=a(i,j,k)
    end do
  end do
end do

kmin=1
if (w_node .and. coord == 0) kmin=2
!$acc parallel loop collapse(3) default(present) async(1)                    &
!$acc private(phix,nxv,nyv,nzv,nmag,x,y,z,dstep,x1,y1,z1,x2,y2,z2,phi1,    &
!$acc         a1,a2)
do k=kmin,nz-1
  do j=1,ny
    do i=1,nx
      if (w_node) then
        phix=0.5_rprec*(phi(i,j,k)+phi(i,j,k-1))
        nxv=0.5_rprec*(normal(1,i,j,k)+normal(1,i,j,k-1))
        nyv=0.5_rprec*(normal(2,i,j,k)+normal(2,i,j,k-1))
        nzv=0.5_rprec*(normal(3,i,j,k)+normal(3,i,j,k-1))
      else
        phix=phi(i,j,k)
        nxv=normal(1,i,j,k)
        nyv=normal(2,i,j,k)
        nzv=normal(3,i,j,k)
      end if
      if (phix >= -phi_cutoff .and. phix < phi_zero) then
        nmag=sqrt(nxv*nxv+nyv*nyv+nzv*nzv)
        if (w_node .and. nmag <= eps) then
          nxv=normal(1,i,j,k)
          nyv=normal(2,i,j,k)
          nzv=normal(3,i,j,k)
          nmag=sqrt(nxv*nxv+nyv*nyv+nzv*nzv)
        end if
        if (nmag > eps) then
          nxv=nxv/nmag
          nyv=nyv/nmag
          nzv=nzv/nmag
          x=real(i-1,rprec)*dx
          y=real(j-1,rprec)*dy
          if (w_node) then
            z=real(k-1,rprec)*dz
          else
            z=(real(k,rprec)-0.5_rprec)*dz
          end if
          dstep=filter_size*dx
          x1=x+dstep*nxv
          y1=y+dstep*nyv
          z1=z+dstep*nzv
#ifdef PPMPI
          if (nproc > 1) then
            phi1=ls_interp_field_halo(phi,phibot,phitop,nphibot,nphitop,   &
                                      lbz,x1,y1,z1,.false.,ld,nx,ny,nz,    &
                                      dx,dy,dz,L_x,L_y)
          else
#endif
            phi1=ls_interp_field(phi,x1,y1,z1,.false.,LS_GRID_ARGS)
#ifdef PPMPI
          end if
#endif
          if (phi1 < phi_zero) then
            dstep=1.5_rprec*dstep
            x1=x+dstep*nxv
            y1=y+dstep*nyv
            z1=z+dstep*nzv
#ifdef PPMPI
            if (nproc > 1) then
              phi1=ls_interp_field_halo(phi,phibot,phitop,nphibot,nphitop,&
                                        lbz,x1,y1,z1,.false.,ld,nx,ny,nz,  &
                                        dx,dy,dz,L_x,L_y)
            else
#endif
              phi1=ls_interp_field(phi,x1,y1,z1,.false.,LS_GRID_ARGS)
#ifdef PPMPI
            end if
#endif
          end if
          if (phi1 < phi_zero) then
            a(i,j,k)=0._rprec
          else
            x2=x1+dstep*nxv
            y2=y1+dstep*nyv
            z2=z1+dstep*nzv
#ifdef PPMPI
            if (nproc > 1) then
              a1=ls_interp_field_halo(tau_source,abot,atop,ntaubot,ntautop,&
                                      lbz,x1,y1,z1,w_node,ld,nx,ny,nz,      &
                                      dx,dy,dz,L_x,L_y)
              a2=ls_interp_field_halo(tau_source,abot,atop,ntaubot,ntautop,&
                                      lbz,x2,y2,z2,w_node,ld,nx,ny,nz,      &
                                      dx,dy,dz,L_x,L_y)
            else
#endif
              a1=ls_interp_field(tau_source,x1,y1,z1,w_node,LS_GRID_ARGS)
              a2=ls_interp_field(tau_source,x2,y2,z2,w_node,LS_GRID_ARGS)
#ifdef PPMPI
            end if
#endif
            a(i,j,k)=2._rprec*a1-a2
          end if
        end if
      end if
    end do
  end do
end do
end subroutine extrap_tau_field_gpu

!*******************************************************************************
subroutine smooth_tau_gpu(phi_limit)
!*******************************************************************************
! Match the CPU multi-color SOR smoother. Each color is independent and runs
! in one full-domain kernel, avoiding both wavefront races and launch-heavy
! per-diagonal synchronization.
real(rprec), intent(in) :: phi_limit
integer, parameter :: niter=5
real(rprec), parameter :: omega=1.5_rprec
integer :: iter,color,ncolors,point_color,cx,cy,i,j,k
integer :: im1,ip1,jm1,jp1,kmin,kmax
real(rprec) :: phi_uv,phi_w,update
logical :: smooth_3d

if (trim(smooth_mode) == 'xy') then
  kmin=1
  kmax=nz
  smooth_3d=.false.
else if (trim(smooth_mode) == '3d') then
  kmin=2
  kmax=nz-1
  smooth_3d=.true.
else
  call error('level_set_gpu_m.smooth_tau_gpu',                             &
      'smooth_mode must be exactly "xy" or "3d": ' // trim(smooth_mode))
end if

ncolors=2
if (modulo(nx,2) /= 0 .or. modulo(ny,2) /= 0) ncolors=3
do iter=1,niter
  do color=0,ncolors-1
    !$acc parallel loop collapse(3) default(present) async(1)               &
    !$acc firstprivate(color,ncolors,smooth_3d,kmin)                        &
    !$acc private(point_color,cx,cy,im1,ip1,jm1,jp1,phi_uv,phi_w,update)
    do k=kmin,kmax
      do j=1,ny
        do i=1,nx
        if (ncolors == 2) then
          point_color=modulo(i+j,2)
          if (smooth_3d) point_color=modulo(point_color+k,2)
        else
          cx=modulo(i-1,2)
          cy=modulo(j-1,2)
          if (modulo(nx,2) /= 0 .and. i == nx) cx=2
          if (modulo(ny,2) /= 0 .and. j == ny) cy=2
          point_color=modulo(cx+cy,3)
          if (smooth_3d) then
            point_color=modulo(point_color+modulo(k-kmin,2),3)
          end if
        end if
        if (point_color == color) then
        im1=i-1
        if (im1 < 1) im1=nx
        ip1=i+1
        if (ip1 > nx) ip1=1
        jm1=j-1
        if (jm1 < 1) jm1=ny
        jp1=j+1
        if (jp1 > ny) jp1=1

        phi_uv=phi(i,j,k)
        if (phi_uv < phi_limit) then
          if (smooth_3d) then
            update=(txx(im1,j,k)+txx(ip1,j,k)+txx(i,jm1,k)+txx(i,jp1,k) + &
                    txx(i,j,k-1)+txx(i,j,k+1))/6._rprec
          else
            update=(txx(im1,j,k)+txx(ip1,j,k)+txx(i,jm1,k)+txx(i,jp1,k)) / &
                    4._rprec
          end if
          txx(i,j,k)=(1._rprec-omega)*txx(i,j,k)+omega*update
          if (smooth_3d) then
            update=(txy(im1,j,k)+txy(ip1,j,k)+txy(i,jm1,k)+txy(i,jp1,k) + &
                    txy(i,j,k-1)+txy(i,j,k+1))/6._rprec
          else
            update=(txy(im1,j,k)+txy(ip1,j,k)+txy(i,jm1,k)+txy(i,jp1,k)) / &
                    4._rprec
          end if
          txy(i,j,k)=(1._rprec-omega)*txy(i,j,k)+omega*update
          if (smooth_3d) then
            update=(tyy(im1,j,k)+tyy(ip1,j,k)+tyy(i,jm1,k)+tyy(i,jp1,k) + &
                    tyy(i,j,k-1)+tyy(i,j,k+1))/6._rprec
          else
            update=(tyy(im1,j,k)+tyy(ip1,j,k)+tyy(i,jm1,k)+tyy(i,jp1,k)) / &
                    4._rprec
          end if
          tyy(i,j,k)=(1._rprec-omega)*tyy(i,j,k)+omega*update
          if (smooth_3d) then
            update=(tzz(im1,j,k)+tzz(ip1,j,k)+tzz(i,jm1,k)+tzz(i,jp1,k) + &
                    tzz(i,j,k-1)+tzz(i,j,k+1))/6._rprec
          else
            update=(tzz(im1,j,k)+tzz(ip1,j,k)+tzz(i,jm1,k)+tzz(i,jp1,k)) / &
                    4._rprec
          end if
          tzz(i,j,k)=(1._rprec-omega)*tzz(i,j,k)+omega*update
        end if

        if (coord == 0 .and. k == 1) then
          phi_w=phi(i,j,k)
        else
          phi_w=0.5_rprec*(phi(i,j,k)+phi(i,j,k-1))
        end if
        if (phi_w < phi_limit) then
          if (smooth_3d) then
            update=(txz(im1,j,k)+txz(ip1,j,k)+txz(i,jm1,k)+txz(i,jp1,k) + &
                    txz(i,j,k-1)+txz(i,j,k+1))/6._rprec
          else
            update=(txz(im1,j,k)+txz(ip1,j,k)+txz(i,jm1,k)+txz(i,jp1,k)) / &
                    4._rprec
          end if
          txz(i,j,k)=(1._rprec-omega)*txz(i,j,k)+omega*update
          if (smooth_3d) then
            update=(tyz(im1,j,k)+tyz(ip1,j,k)+tyz(i,jm1,k)+tyz(i,jp1,k) + &
                    tyz(i,j,k-1)+tyz(i,j,k+1))/6._rprec
          else
            update=(tyz(im1,j,k)+tyz(ip1,j,k)+tyz(i,jm1,k)+tyz(i,jp1,k)) / &
                    4._rprec
          end if
          tyz(i,j,k)=(1._rprec-omega)*tyz(i,j,k)+omega*update
        end if
        end if
        end do
      end do
    end do
  end do
end do
end subroutine smooth_tau_gpu

!*******************************************************************************
subroutine smooth_field_gpu(a, albz, w_node, phi_limit)
!*******************************************************************************
integer, intent(in) :: albz
real(rprec), intent(inout) :: a(ld,ny,albz:nz)
logical, intent(in) :: w_node
real(rprec), intent(in) :: phi_limit
integer, parameter :: niter=5
real(rprec), parameter :: omega=1.5_rprec
integer :: iter,color,ncolors,point_color,cx,cy,i,j,k
integer :: im1,ip1,jm1,jp1,kmin,kmax,shift
real(rprec) :: phiv,update
logical :: smooth_3d

shift=0
if (w_node) shift=1

if (trim(smooth_mode) == 'xy') then
  kmin=1
  kmax=nz
  smooth_3d=.false.
else if (trim(smooth_mode) == '3d') then
  kmin=2
  kmax=nz-1
  smooth_3d=.true.
else
  call error('level_set_gpu_m.smooth_field_gpu',                           &
      'smooth_mode must be exactly "xy" or "3d": ' // trim(smooth_mode))
end if

ncolors=2
if (modulo(nx,2) /= 0 .or. modulo(ny,2) /= 0) ncolors=3
do iter=1,niter
  do color=0,ncolors-1
    !$acc parallel loop collapse(3) default(present) async(1)               &
    !$acc firstprivate(color,ncolors,smooth_3d,kmin)                         &
    !$acc private(point_color,cx,cy,im1,ip1,jm1,jp1,phiv,update)
    do k=kmin,kmax
      do j=1,ny
        do i=1,nx
          if (ncolors == 2) then
            point_color=modulo(i+j,2)
            if (smooth_3d) point_color=modulo(point_color+k,2)
          else
            cx=modulo(i-1,2)
            cy=modulo(j-1,2)
            if (modulo(nx,2) /= 0 .and. i == nx) cx=2
            if (modulo(ny,2) /= 0 .and. j == ny) cy=2
            point_color=modulo(cx+cy,3)
            if (smooth_3d) then
              point_color=modulo(point_color+modulo(k-kmin,2),3)
            end if
          end if
          if (point_color == color) then
          if (coord == 0 .and. k == shift) then
            phiv=phi(i,j,k)
          else
            phiv=0.5_rprec*(phi(i,j,k)+phi(i,j,k-shift))
          end if
          if (phiv < phi_limit) then
            im1=i-1
            if (im1 < 1) im1=nx
            ip1=i+1
            if (ip1 > nx) ip1=1
            jm1=j-1
            if (jm1 < 1) jm1=ny
            jp1=j+1
            if (jp1 > ny) jp1=1
            if (smooth_3d) then
              update=(a(im1,j,k)+a(ip1,j,k)+a(i,jm1,k)+a(i,jp1,k) +       &
                      a(i,j,k-1)+a(i,j,k+1))/6._rprec
            else
              update=(a(im1,j,k)+a(ip1,j,k)+a(i,jm1,k)+a(i,jp1,k)) /     &
                      4._rprec
            end if
            a(i,j,k)=(1._rprec-omega)*a(i,j,k)+omega*update
          end if
          end if
        end do
      end do
    end do
  end do
end do
end subroutine smooth_field_gpu

!*******************************************************************************
subroutine smooth_lag_fields_gpu(S11,S12,S13,S22,S23,S33)
!*******************************************************************************
! Batch the nine independent fields used by model 4/5 into one color sequence.
! This preserves each field's multi-color SOR order while avoiding nine copies
! of the kernel-launch sequence.
real(rprec), intent(inout) :: S11(ld,ny,nz), S12(ld,ny,nz), S13(ld,ny,nz)
real(rprec), intent(inout) :: S22(ld,ny,nz), S23(ld,ny,nz), S33(ld,ny,nz)
integer, parameter :: niter=5
real(rprec), parameter :: omega=1.5_rprec
integer :: iter,color,ncolors,point_color,cx,cy,i,j,k
integer :: im1,ip1,jm1,jp1,kmin,kmax
real(rprec) :: phi_uv,phi_w,update
logical :: smooth_3d

if (trim(smooth_mode) == 'xy') then
  kmin=1
  kmax=nz
  smooth_3d=.false.
else if (trim(smooth_mode) == '3d') then
  kmin=2
  kmax=nz-1
  smooth_3d=.true.
else
  call error('level_set_gpu_m.smooth_lag_fields_gpu',                      &
      'smooth_mode must be exactly "xy" or "3d": ' // trim(smooth_mode))
end if

ncolors=2
if (modulo(nx,2) /= 0 .or. modulo(ny,2) /= 0) ncolors=3
do iter=1,niter
  do color=0,ncolors-1
    !$acc parallel loop collapse(3) default(present) async(1)               &
    !$acc firstprivate(color,ncolors,smooth_3d,kmin)                        &
    !$acc private(point_color,cx,cy,im1,ip1,jm1,jp1,phi_uv,phi_w,update)
    do k=kmin,kmax
      do j=1,ny
        do i=1,nx
          if (ncolors == 2) then
            point_color=modulo(i+j,2)
            if (smooth_3d) point_color=modulo(point_color+k,2)
          else
            cx=modulo(i-1,2)
            cy=modulo(j-1,2)
            if (modulo(nx,2) /= 0 .and. i == nx) cx=2
            if (modulo(ny,2) /= 0 .and. j == ny) cy=2
            point_color=modulo(cx+cy,3)
            if (smooth_3d) then
              point_color=modulo(point_color+modulo(k-kmin,2),3)
            end if
          end if
          if (point_color /= color) cycle

          im1=i-1
          if (im1 < 1) im1=nx
          ip1=i+1
          if (ip1 > nx) ip1=1
          jm1=j-1
          if (jm1 < 1) jm1=ny
          jp1=j+1
          if (jp1 > ny) jp1=1

          phi_uv=phi(i,j,k)
          if (phi_uv < 0._rprec) then
            if (smooth_3d) then
              update=(u(im1,j,k)+u(ip1,j,k)+u(i,jm1,k)+u(i,jp1,k) +       &
                      u(i,j,k-1)+u(i,j,k+1))/6._rprec
            else
              update=(u(im1,j,k)+u(ip1,j,k)+u(i,jm1,k)+u(i,jp1,k)) /     &
                      4._rprec
            end if
            u(i,j,k)=(1._rprec-omega)*u(i,j,k)+omega*update
            if (smooth_3d) then
              update=(v(im1,j,k)+v(ip1,j,k)+v(i,jm1,k)+v(i,jp1,k) +       &
                      v(i,j,k-1)+v(i,j,k+1))/6._rprec
            else
              update=(v(im1,j,k)+v(ip1,j,k)+v(i,jm1,k)+v(i,jp1,k)) /     &
                      4._rprec
            end if
            v(i,j,k)=(1._rprec-omega)*v(i,j,k)+omega*update
          end if

          if (coord == 0 .and. k == 1) then
            phi_w=phi(i,j,k)
          else
            phi_w=0.5_rprec*(phi(i,j,k)+phi(i,j,k-1))
          end if
          if (phi_w < 0._rprec) then
            if (smooth_3d) then
              update=(w(im1,j,k)+w(ip1,j,k)+w(i,jm1,k)+w(i,jp1,k) +       &
                      w(i,j,k-1)+w(i,j,k+1))/6._rprec
            else
              update=(w(im1,j,k)+w(ip1,j,k)+w(i,jm1,k)+w(i,jp1,k)) /     &
                      4._rprec
            end if
            w(i,j,k)=(1._rprec-omega)*w(i,j,k)+omega*update
            if (smooth_3d) then
              update=(S11(im1,j,k)+S11(ip1,j,k)+S11(i,jm1,k)+            &
                      S11(i,jp1,k)+S11(i,j,k-1)+S11(i,j,k+1))/6._rprec
            else
              update=(S11(im1,j,k)+S11(ip1,j,k)+S11(i,jm1,k)+            &
                      S11(i,jp1,k))/4._rprec
            end if
            S11(i,j,k)=(1._rprec-omega)*S11(i,j,k)+omega*update
            if (smooth_3d) then
              update=(S12(im1,j,k)+S12(ip1,j,k)+S12(i,jm1,k)+            &
                      S12(i,jp1,k)+S12(i,j,k-1)+S12(i,j,k+1))/6._rprec
            else
              update=(S12(im1,j,k)+S12(ip1,j,k)+S12(i,jm1,k)+            &
                      S12(i,jp1,k))/4._rprec
            end if
            S12(i,j,k)=(1._rprec-omega)*S12(i,j,k)+omega*update
            if (smooth_3d) then
              update=(S13(im1,j,k)+S13(ip1,j,k)+S13(i,jm1,k)+            &
                      S13(i,jp1,k)+S13(i,j,k-1)+S13(i,j,k+1))/6._rprec
            else
              update=(S13(im1,j,k)+S13(ip1,j,k)+S13(i,jm1,k)+            &
                      S13(i,jp1,k))/4._rprec
            end if
            S13(i,j,k)=(1._rprec-omega)*S13(i,j,k)+omega*update
            if (smooth_3d) then
              update=(S22(im1,j,k)+S22(ip1,j,k)+S22(i,jm1,k)+            &
                      S22(i,jp1,k)+S22(i,j,k-1)+S22(i,j,k+1))/6._rprec
            else
              update=(S22(im1,j,k)+S22(ip1,j,k)+S22(i,jm1,k)+            &
                      S22(i,jp1,k))/4._rprec
            end if
            S22(i,j,k)=(1._rprec-omega)*S22(i,j,k)+omega*update
            if (smooth_3d) then
              update=(S23(im1,j,k)+S23(ip1,j,k)+S23(i,jm1,k)+            &
                      S23(i,jp1,k)+S23(i,j,k-1)+S23(i,j,k+1))/6._rprec
            else
              update=(S23(im1,j,k)+S23(ip1,j,k)+S23(i,jm1,k)+            &
                      S23(i,jp1,k))/4._rprec
            end if
            S23(i,j,k)=(1._rprec-omega)*S23(i,j,k)+omega*update
            if (smooth_3d) then
              update=(S33(im1,j,k)+S33(ip1,j,k)+S33(i,jm1,k)+            &
                      S33(i,jp1,k)+S33(i,j,k-1)+S33(i,j,k+1))/6._rprec
            else
              update=(S33(im1,j,k)+S33(ip1,j,k)+S33(i,jm1,k)+            &
                      S33(i,jp1,k))/4._rprec
            end if
            S33(i,j,k)=(1._rprec-omega)*S33(i,j,k)+omega*update
          end if
        end do
      end do
    end do
  end do
end do
end subroutine smooth_lag_fields_gpu

!*******************************************************************************
subroutine level_set_lag_dyn_gpu_core(S11,S12,S13,S22,S23,S33,normal,      &
                                      modify_beta_enabled)
!*******************************************************************************
! Device implementation of level_set_lag_dyn. This is the Level Set
! preconditioning required by the Lagrangian similarity model (sgs_model=4),
! before its F-history interpolation and coefficient update.
use param, only : lbc_mom
use sgs_param, only : F_LM, F_MM, beta
use test_filtermodule, only : filter_size
real(rprec), intent(inout) :: S11(ld,ny,nz), S12(ld,ny,nz), S13(ld,ny,nz), &
                              S22(ld,ny,nz), S23(ld,ny,nz), S33(ld,ny,nz)
real(rprec), intent(in) :: normal(3,ld,ny,lbz:nz)
logical, intent(in) :: modify_beta_enabled
integer :: i,j,k,s
real(rprec), parameter :: c1=0.65_rprec, c2=0.7_rprec
real(rprec) :: phix,phi1,dphi,x,y,z,xp,yp,zp,nxv,nyv,nzv
real(rprec) :: delta_filter,dmin,zglobal

call smooth_lag_fields_gpu(S11,S12,S13,S22,S23,S33)

! F_LM is inactive inside the solid and in its filter-width buffer.
!$acc parallel loop collapse(3) default(present) async(1) private(s,phix)
do k=1,nz-1
  do j=1,ny
    do i=1,nx
      s=1
      if (coord == 0 .and. k == 1) s=0
      phix=0.5_rprec*(phi(i,j,k)+phi(i,j,k-s))
      if (phix < filter_size*dx) F_LM(i,j,k)=0._rprec
    end do
  end do
end do

! Freeze F_MM so every near-interface target sees the same source field.
!$acc parallel loop collapse(3) present(F_MM,fmm_source) async(1)
do k=1,nz
  do j=1,ny
    do i=1,ld
      fmm_source(i,j,k)=F_MM(i,j,k)
    end do
  end do
end do

dphi=sqrt(dx*dx+dy*dy+dz*dz)
phi1=-dphi
!$acc parallel loop collapse(3) default(present) async(1)                    &
!$acc private(s,phix,x,y,z,nxv,nyv,nzv,xp,yp,zp)
do k=1,nz-1
  do j=1,ny
    do i=1,nx
      s=1
      if (coord == 0 .and. k == 1) s=0
      phix=0.5_rprec*(phi(i,j,k)+phi(i,j,k-s))
      if (phi1 < phix .and. phix < 0._rprec) then
        x=real(i-1,rprec)*dx
        y=real(j-1,rprec)*dy
        z=(real(k,rprec)-0.5_rprec*real(1+s,rprec))*dz
        nxv=0.5_rprec*(normal(1,i,j,k)+normal(1,i,j,k-s))
        nyv=0.5_rprec*(normal(2,i,j,k)+normal(2,i,j,k-s))
        nzv=0.5_rprec*(normal(3,i,j,k)+normal(3,i,j,k-s))
        xp=x+dphi*nxv
        yp=y+dphi*nyv
        zp=z+dphi*nzv
#ifdef PPMPI
        if (nproc > 1) then
          F_MM(i,j,k)=ls_interp_field_halo(fmm_source,FMMbot,FMMtop,        &
                                           nFMMbot,nFMMtop,1,xp,yp,zp,      &
                                           .true.,ld,nx,ny,nz,              &
                                           dx,dy,dz,L_x,L_y)
        else
#endif
          F_MM(i,j,k)=ls_interp_field(fmm_source,xp,yp,zp,.true.,           &
                                      ld,nx,ny,nz,1,dx,dy,dz,L_x,L_y)
#ifdef PPMPI
        end if
#endif
      end if
    end do
  end do
end do

if (modify_beta_enabled) then
  ! Match modify_beta(): retain beta outside the solid and account for the
  ! physical lower wall when it is active.
  delta_filter=filter_size*(dx*dy*dz)**(1._rprec/3._rprec)
  !$acc parallel loop collapse(3) default(present) async(1)                  &
  !$acc private(s,phix,zglobal,dmin)
  do k=1,nz
    do j=1,ny
      do i=1,nx
        s=1
        if (coord == 0 .and. k == 1) s=0
        phix=0.5_rprec*(phi(i,j,k)+phi(i,j,k-s))
        if (phix > 0._rprec) then
          if (lbc_mom == 0) then
            dmin=phix
          else
            zglobal=real(coord*(nz-1)+k-1,rprec)*dz
            dmin=min(zglobal,phix)
          end if
          beta(i,j,k)=1._rprec-c1*exp(-c2*dmin/delta_filter)
        end if
      end do
    end do
  end do
end if
end subroutine level_set_lag_dyn_gpu_core

end module level_set_gpu_m
#endif
