!*******************************************************************************
module hit_inflow_gpu
!*******************************************************************************
! GPU helper module for optional HIT inflow support.
!
! Ownership map:
!   - hit_gpu_setup: copy HIT field/table data into GPU-resident storage
!   - hit_fringe_setup_gpu: copy fringe coefficients for GPU use
!   - hit_compute_plane_gpu: interpolate one inflow plane on the GPU
!   - hit_apply_fringe_gpu: apply GPU fringe blending
!
! This module is compiled only for the optional `USE_HIT` path.  It is not part
! of the validated wind-farm production configuration unless HIT inflow is
! explicitly enabled and validated.
use types, only : rprec
use param, only : ld, ny, nz, lbz, dy, dz
#ifdef PPMPI
use param, only : coord
#endif
#ifdef ENABLE_CUDA
use cudafor
#endif
implicit none

private
public :: hit_cuda_enabled, hit_gpu_setup, hit_fringe_setup_gpu,              &
    hit_compute_plane_gpu, hit_apply_fringe_gpu

#ifdef ENABLE_CUDA
real(rprec), managed, save, allocatable, dimension(:) :: hit_x_cuda,          &
    hit_y_cuda, hit_z_cuda
real(rprec), managed, save, allocatable, dimension(:,:,:) :: hit_u_cuda,      &
    hit_v_cuda, hit_w_cuda
real(rprec), managed, save, allocatable, dimension(:,:) :: hit_u_plane_cuda,  &
    hit_v_plane_cuda, hit_w_plane_cuda
integer, managed, save, allocatable, dimension(:) :: hit_iwrap_cuda
real(rprec), managed, save, allocatable, dimension(:) :: hit_alpha_cuda,      &
    hit_beta_cuda
#elif defined(PPLES_GPU)
real(rprec), save, allocatable, dimension(:) :: hit_x_cuda, hit_y_cuda,       &
    hit_z_cuda
real(rprec), save, allocatable, dimension(:,:,:) :: hit_u_cuda, hit_v_cuda,   &
    hit_w_cuda
real(rprec), save, allocatable, dimension(:,:) :: hit_u_plane_cuda,           &
    hit_v_plane_cuda, hit_w_plane_cuda
integer, save, allocatable, dimension(:) :: hit_iwrap_cuda
real(rprec), save, allocatable, dimension(:) :: hit_alpha_cuda, hit_beta_cuda
!$acc declare create(hit_x_cuda, hit_y_cuda, hit_z_cuda)
!$acc declare create(hit_u_cuda, hit_v_cuda, hit_w_cuda)
!$acc declare create(hit_u_plane_cuda, hit_v_plane_cuda, hit_w_plane_cuda)
!$acc declare create(hit_iwrap_cuda, hit_alpha_cuda, hit_beta_cuda)
#endif

#if defined(ENABLE_CUDA) || defined(PPLES_GPU)
integer, save :: nx_hit_cuda = 0, ny_hit_cuda = 0, nz_hit_cuda = 0
integer, save :: hit_fringe_nx_cuda = 0
#endif

contains

!*******************************************************************************
logical function hit_cuda_enabled()
!*******************************************************************************
implicit none

#if defined(ENABLE_CUDA) || defined(PPLES_GPU)
hit_cuda_enabled = .true.
#else
hit_cuda_enabled = .false.
#endif

end function hit_cuda_enabled

!*******************************************************************************
subroutine hit_gpu_setup(nx_hit, ny_hit, nz_hit, x, y, z, u, v, w)
!*******************************************************************************
integer, intent(in) :: nx_hit, ny_hit, nz_hit
real(rprec), intent(in) :: x(nx_hit), y(ny_hit), z(nz_hit)
real(rprec), intent(in) :: u(nx_hit, ny_hit, nz_hit),                         &
    v(nx_hit, ny_hit, nz_hit), w(nx_hit, ny_hit, nz_hit)

#if defined(ENABLE_CUDA) || defined(PPLES_GPU)
if (allocated(hit_x_cuda)) deallocate(hit_x_cuda)
if (allocated(hit_y_cuda)) deallocate(hit_y_cuda)
if (allocated(hit_z_cuda)) deallocate(hit_z_cuda)
if (allocated(hit_u_cuda)) deallocate(hit_u_cuda)
if (allocated(hit_v_cuda)) deallocate(hit_v_cuda)
if (allocated(hit_w_cuda)) deallocate(hit_w_cuda)
if (allocated(hit_u_plane_cuda)) deallocate(hit_u_plane_cuda)
if (allocated(hit_v_plane_cuda)) deallocate(hit_v_plane_cuda)
if (allocated(hit_w_plane_cuda)) deallocate(hit_w_plane_cuda)

nx_hit_cuda = nx_hit
ny_hit_cuda = ny_hit
nz_hit_cuda = nz_hit

allocate(hit_x_cuda(nx_hit))
allocate(hit_y_cuda(ny_hit))
allocate(hit_z_cuda(nz_hit))
allocate(hit_u_cuda(nx_hit, ny_hit, nz_hit))
allocate(hit_v_cuda(nx_hit, ny_hit, nz_hit))
allocate(hit_w_cuda(nx_hit, ny_hit, nz_hit))
allocate(hit_u_plane_cuda(ny, nz))
allocate(hit_v_plane_cuda(ny, nz))
allocate(hit_w_plane_cuda(ny, nz))

hit_x_cuda = x
hit_y_cuda = y
hit_z_cuda = z
hit_u_cuda = u
hit_v_cuda = v
hit_w_cuda = w
#endif

#ifdef PPLES_GPU
!$acc update device(hit_x_cuda, hit_y_cuda, hit_z_cuda)
!$acc update device(hit_u_cuda, hit_v_cuda, hit_w_cuda)
#endif

end subroutine hit_gpu_setup

!*******************************************************************************
subroutine hit_fringe_setup_gpu(fringe_nx, iwrap, alpha, beta)
!*******************************************************************************
integer, intent(in) :: fringe_nx
integer, intent(in) :: iwrap(fringe_nx)
real(rprec), intent(in) :: alpha(fringe_nx), beta(fringe_nx)

#if defined(ENABLE_CUDA) || defined(PPLES_GPU)
if (allocated(hit_iwrap_cuda)) deallocate(hit_iwrap_cuda)
if (allocated(hit_alpha_cuda)) deallocate(hit_alpha_cuda)
if (allocated(hit_beta_cuda)) deallocate(hit_beta_cuda)

hit_fringe_nx_cuda = fringe_nx
allocate(hit_iwrap_cuda(fringe_nx))
allocate(hit_alpha_cuda(fringe_nx))
allocate(hit_beta_cuda(fringe_nx))
hit_iwrap_cuda = iwrap
hit_alpha_cuda = alpha
hit_beta_cuda = beta
#endif

#ifdef PPLES_GPU
!$acc update device(hit_iwrap_cuda, hit_alpha_cuda, hit_beta_cuda)
#endif

end subroutine hit_fringe_setup_gpu

!*******************************************************************************
subroutine hit_compute_plane_gpu(xloc, ti_out, inflow_velocity)
!*******************************************************************************
real(rprec), intent(in) :: xloc, ti_out, inflow_velocity

#if defined(ENABLE_CUDA) || defined(PPLES_GPU)
integer :: i0, i1, j0, j1, k0, k1
integer :: i0_x, i1_x
integer :: i, j, k
real(rprec) :: xp, yp, zp, zwp
real(rprec) :: xd, xd_x, yd, zd
real(rprec) :: x0, x1, y0, y1, z0, z1
real(rprec) :: c00, c01, c10, c11, c0, c1
real(rprec) :: inv_dx_hit, inv_dy_hit, inv_dz_hit

if (.not. allocated(hit_x_cuda)) return

if (nx_hit_cuda > 1) then
    inv_dx_hit = 1._rprec / (hit_x_cuda(nx_hit_cuda) - hit_x_cuda(1))           &
        * real(nx_hit_cuda - 1, rprec)
else
    inv_dx_hit = 0._rprec
endif
if (ny_hit_cuda > 1) then
    inv_dy_hit = 1._rprec / (hit_y_cuda(ny_hit_cuda) - hit_y_cuda(1))           &
        * real(ny_hit_cuda - 1, rprec)
else
    inv_dy_hit = 0._rprec
endif
if (nz_hit_cuda > 1) then
    inv_dz_hit = 1._rprec / (hit_z_cuda(nz_hit_cuda) - hit_z_cuda(1))           &
        * real(nz_hit_cuda - 1, rprec)
else
    inv_dz_hit = 0._rprec
endif

xp = xloc
if (xp <= hit_x_cuda(1) .or. nx_hit_cuda <= 1) then
    i0_x = 1
    i1_x = 1
    xd_x = 1._rprec
else if (xp >= hit_x_cuda(nx_hit_cuda)) then
    i0_x = nx_hit_cuda
    i1_x = nx_hit_cuda
    xd_x = 1._rprec
else
    i0_x = max(1, min(nx_hit_cuda - 1,                                      &
        int((xp - hit_x_cuda(1)) * inv_dx_hit) + 1))
    i1_x = i0_x + 1
    x0 = hit_x_cuda(i0_x)
    x1 = hit_x_cuda(i1_x)
    xd_x = (xp - x0) / (x1 - x0)
endif
#endif

#ifdef PPLES_GPU
!$acc parallel loop collapse(2) default(present) private(i0,i1,j0,j1,k0,k1,    &
!$acc& xp,yp,zp,zwp,xd,yd,zd,x0,x1,y0,y1,z0,z1,c00,c01,c10,c11,c0,c1)
do k = 1, nz
    do j = 1, ny
#elif defined(ENABLE_CUDA)
!$cuf kernel do(2) <<<*, *>>>
do k = 1, nz
    do j = 1, ny
#endif
#if defined(ENABLE_CUDA) || defined(PPLES_GPU)
        i0 = i0_x
        i1 = i1_x
        xd = xd_x
        yp = (j - 1) * dy
#ifdef PPMPI
        zp = (coord * (nz - 1) + k - 0.5_rprec) * dz
#else
        zp = (k - 0.5_rprec) * dz
#endif
        zwp = zp - 0.5_rprec * dz

        if (yp <= hit_y_cuda(1) .or. ny_hit_cuda <= 1) then
            j0 = 1
            j1 = 1
            yd = 1._rprec
        else if (yp >= hit_y_cuda(ny_hit_cuda)) then
            j0 = ny_hit_cuda
            j1 = ny_hit_cuda
            yd = 1._rprec
        else
            j0 = max(1, min(ny_hit_cuda - 1,                                 &
                int((yp - hit_y_cuda(1)) * inv_dy_hit) + 1))
            j1 = j0 + 1
            y0 = hit_y_cuda(j0)
            y1 = hit_y_cuda(j1)
            yd = (yp - y0) / (y1 - y0)
        endif

        if (zp <= hit_z_cuda(1) .or. nz_hit_cuda <= 1) then
            k0 = 1
            k1 = 1
            zd = 1._rprec
        else if (zp >= hit_z_cuda(nz_hit_cuda)) then
            k0 = nz_hit_cuda
            k1 = nz_hit_cuda
            zd = 1._rprec
        else
            k0 = max(1, min(nz_hit_cuda - 1,                                 &
                int((zp - hit_z_cuda(1)) * inv_dz_hit) + 1))
            k1 = k0 + 1
            z0 = hit_z_cuda(k0)
            z1 = hit_z_cuda(k1)
            zd = (zp - z0) / (z1 - z0)
        endif

        c00 = hit_u_cuda(i0,j0,k0) * (1._rprec - xd) +                       &
            hit_u_cuda(i1,j0,k0) * xd
        c01 = hit_u_cuda(i0,j0,k1) * (1._rprec - xd) +                       &
            hit_u_cuda(i1,j0,k1) * xd
        c10 = hit_u_cuda(i0,j1,k0) * (1._rprec - xd) +                       &
            hit_u_cuda(i1,j1,k0) * xd
        c11 = hit_u_cuda(i0,j1,k1) * (1._rprec - xd) +                       &
            hit_u_cuda(i1,j1,k1) * xd
        c0 = c00 * (1._rprec - yd) + c10 * yd
        c1 = c01 * (1._rprec - yd) + c11 * yd
        hit_u_plane_cuda(j,k) = inflow_velocity * (1._rprec + ti_out *        &
            (c0 * (1._rprec - zd) + c1 * zd))

        c00 = hit_v_cuda(i0,j0,k0) * (1._rprec - xd) +                       &
            hit_v_cuda(i1,j0,k0) * xd
        c01 = hit_v_cuda(i0,j0,k1) * (1._rprec - xd) +                       &
            hit_v_cuda(i1,j0,k1) * xd
        c10 = hit_v_cuda(i0,j1,k0) * (1._rprec - xd) +                       &
            hit_v_cuda(i1,j1,k0) * xd
        c11 = hit_v_cuda(i0,j1,k1) * (1._rprec - xd) +                       &
            hit_v_cuda(i1,j1,k1) * xd
        c0 = c00 * (1._rprec - yd) + c10 * yd
        c1 = c01 * (1._rprec - yd) + c11 * yd
        hit_v_plane_cuda(j,k) = inflow_velocity * ti_out *                   &
            (c0 * (1._rprec - zd) + c1 * zd)

        if (zwp <= hit_z_cuda(1) .or. nz_hit_cuda <= 1) then
            k0 = 1
            k1 = 1
            zd = 1._rprec
        else if (zwp >= hit_z_cuda(nz_hit_cuda)) then
            k0 = nz_hit_cuda
            k1 = nz_hit_cuda
            zd = 1._rprec
        else
            k0 = max(1, min(nz_hit_cuda - 1,                                 &
                int((zwp - hit_z_cuda(1)) * inv_dz_hit) + 1))
            k1 = k0 + 1
            z0 = hit_z_cuda(k0)
            z1 = hit_z_cuda(k1)
            zd = (zwp - z0) / (z1 - z0)
        endif

        c00 = hit_w_cuda(i0,j0,k0) * (1._rprec - xd) +                       &
            hit_w_cuda(i1,j0,k0) * xd
        c01 = hit_w_cuda(i0,j0,k1) * (1._rprec - xd) +                       &
            hit_w_cuda(i1,j0,k1) * xd
        c10 = hit_w_cuda(i0,j1,k0) * (1._rprec - xd) +                       &
            hit_w_cuda(i1,j1,k0) * xd
        c11 = hit_w_cuda(i0,j1,k1) * (1._rprec - xd) +                       &
            hit_w_cuda(i1,j1,k1) * xd
        c0 = c00 * (1._rprec - yd) + c10 * yd
        c1 = c01 * (1._rprec - yd) + c11 * yd
        hit_w_plane_cuda(j,k) = inflow_velocity * ti_out *                   &
            (c0 * (1._rprec - zd) + c1 * zd)
    enddo
enddo
#endif

end subroutine hit_compute_plane_gpu

!*******************************************************************************
subroutine hit_apply_fringe_gpu(fringe_nx, iwrap, alpha, beta, u, v, w)
!*******************************************************************************
integer, intent(in) :: fringe_nx
integer, intent(in) :: iwrap(fringe_nx)
real(rprec), intent(in) :: alpha(fringe_nx), beta(fringe_nx)
real(rprec), intent(inout) :: u(ld,ny,lbz:nz), v(ld,ny,lbz:nz), w(ld,ny,lbz:nz)

#if defined(ENABLE_CUDA) || defined(PPLES_GPU)
integer :: i, j, k, i_w
real(rprec) :: alpha_i, beta_i

if (.not. allocated(hit_u_plane_cuda)) return
#endif

#ifdef PPLES_GPU
if (.not. allocated(hit_iwrap_cuda)) return
!$acc parallel loop collapse(3) default(present) private(i_w, alpha_i, beta_i)
do k = 1, nz
    do j = 1, ny
        do i = 1, hit_fringe_nx_cuda
            i_w = hit_iwrap_cuda(i)
            alpha_i = hit_alpha_cuda(i)
            beta_i = hit_beta_cuda(i)
            u(i_w,j,k) = alpha_i * u(i_w,j,k) + beta_i * hit_u_plane_cuda(j,k)
            v(i_w,j,k) = alpha_i * v(i_w,j,k) + beta_i * hit_v_plane_cuda(j,k)
            w(i_w,j,k) = alpha_i * w(i_w,j,k) + beta_i * hit_w_plane_cuda(j,k)
        enddo
    enddo
enddo
#elif defined(ENABLE_CUDA)
attributes(device) :: iwrap, alpha, beta, u, v, w

!$cuf kernel do(3) <<<*, *>>>
do k = 1, nz
    do j = 1, ny
        do i = 1, fringe_nx
            i_w = iwrap(i)
            alpha_i = alpha(i)
            beta_i = beta(i)
            u(i_w,j,k) = alpha_i * u(i_w,j,k) + beta_i * hit_u_plane_cuda(j,k)
            v(i_w,j,k) = alpha_i * v(i_w,j,k) + beta_i * hit_v_plane_cuda(j,k)
            w(i_w,j,k) = alpha_i * w(i_w,j,k) + beta_i * hit_w_plane_cuda(j,k)
        enddo
    enddo
enddo
#endif

end subroutine hit_apply_fringe_gpu

end module hit_inflow_gpu
