program validate_filt_da_cufft
use cudafor
use cufft
implicit none
include 'fftw3.f'

integer, parameter :: dp = kind(1.0d0)
integer, parameter :: nx = 32
integer, parameter :: ny = 24
integer, parameter :: lh = nx / 2 + 1
integer, parameter :: ld = 2 * lh
integer, parameter :: nbatch = 3
real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
real(dp), parameter :: tol = 1.0e-10_dp

integer*8 :: forw, back
real(dp) :: field(ld, ny, nbatch)
real(dp) :: f_cpu(ld, ny, nbatch)
real(dp) :: dfdx_cpu(ld, ny, nbatch)
real(dp) :: dfdy_cpu(ld, ny, nbatch)
real(dp) :: f_gpu(ld, ny, nbatch)
real(dp) :: dfdx_gpu(ld, ny, nbatch)
real(dp) :: dfdy_gpu(ld, ny, nbatch)
real(dp) :: kx(lh, ny)
real(dp) :: ky(lh, ny)
real(dp) :: worst

call init_field(field)
call init_wavenumbers(kx, ky)

call dfftw_plan_dft_r2c_2d(forw, nx, ny, f_cpu(:,:,1), f_cpu(:,:,1),          &
    FFTW_ESTIMATE, FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_2d(back, nx, ny, f_cpu(:,:,1), f_cpu(:,:,1),          &
    FFTW_ESTIMATE, FFTW_UNALIGNED)

call filt_da_cpu(field, f_cpu, dfdx_cpu, dfdy_cpu, kx, ky, forw, back)
call filt_da_gpu_planmany(field, f_gpu, dfdx_gpu, dfdy_gpu, kx, ky)

worst = 0.0_dp
call check_array('filtered field', f_cpu, f_gpu, worst)
call check_array('dfdx', dfdx_cpu, dfdx_gpu, worst)
call check_array('dfdy', dfdy_cpu, dfdy_gpu, worst)

call dfftw_destroy_plan(forw)
call dfftw_destroy_plan(back)

if (worst > tol) then
    write(*,'(a,es12.4,a,es12.4)') 'FAIL: max difference ', worst,            &
        ' exceeds tolerance ', tol
    stop 1
end if

write(*,'(a,es12.4)') 'PASS: FFTW and cuFFT filt_da paths agree; worst diff=', &
    worst

contains

subroutine init_field(a)
implicit none
real(dp), intent(out) :: a(ld, ny, nbatch)
integer :: i, j, k
real(dp) :: x, y

a = 0.0_dp
do k = 1, nbatch
do j = 1, ny
do i = 1, nx
    x = real(i - 1, dp) / real(nx, dp)
    y = real(j - 1, dp) / real(ny, dp)
    a(i,j,k) = sin(2.0_dp*pi*(real(k, dp) + 1.0_dp)*x)                       &
        + 0.5_dp*cos(2.0_dp*pi*3.0_dp*y)                                      &
        + 0.25_dp*sin(2.0_dp*pi*(2.0_dp*x - real(k + 1, dp)*y))
end do
end do
end do

end subroutine init_field

subroutine init_wavenumbers(kx_out, ky_out)
implicit none
real(dp), intent(out) :: kx_out(lh, ny)
real(dp), intent(out) :: ky_out(lh, ny)
integer :: i, j, jm

do j = 1, ny
    if (j <= ny/2 + 1) then
        jm = j - 1
    else
        jm = j - ny - 1
    end if

    do i = 1, lh
        kx_out(i,j) = 2.0_dp*pi*real(i - 1, dp)
        ky_out(i,j) = 2.0_dp*pi*real(jm, dp)
    end do
end do

end subroutine init_wavenumbers

subroutine filt_da_cpu(field_in, f_out, dfdx_out, dfdy_out, kx_in, ky_in,      &
    forw_plan, back_plan)
implicit none
real(dp), intent(in) :: field_in(ld, ny, nbatch)
real(dp), intent(out) :: f_out(ld, ny, nbatch)
real(dp), intent(out) :: dfdx_out(ld, ny, nbatch)
real(dp), intent(out) :: dfdy_out(ld, ny, nbatch)
real(dp), intent(in) :: kx_in(lh, ny)
real(dp), intent(in) :: ky_in(lh, ny)
integer*8, intent(in) :: forw_plan, back_plan
integer :: i, j, k, ir, ii
real(dp) :: const

const = 1.0_dp / real(nx * ny, dp)
f_out = 0.0_dp
dfdx_out = 0.0_dp
dfdy_out = 0.0_dp

do k = 1, nbatch
    f_out(:,:,k) = const * field_in(:,:,k)
    call dfftw_execute_dft_r2c(forw_plan, f_out(:,:,k), f_out(:,:,k))

    f_out(ld-1:ld,:,k) = 0.0_dp
    f_out(:,ny/2+1,k) = 0.0_dp

    do j = 1, ny
    do i = 1, lh
        ir = 2 * i - 1
        ii = 2 * i
        dfdx_out(ir,j,k) = -f_out(ii,j,k) * kx_in(i,j)
        dfdx_out(ii,j,k) =  f_out(ir,j,k) * kx_in(i,j)
        dfdy_out(ir,j,k) = -f_out(ii,j,k) * ky_in(i,j)
        dfdy_out(ii,j,k) =  f_out(ir,j,k) * ky_in(i,j)
    end do
    end do

    call dfftw_execute_dft_c2r(back_plan, f_out(:,:,k), f_out(:,:,k))
    call dfftw_execute_dft_c2r(back_plan, dfdx_out(:,:,k), dfdx_out(:,:,k))
    call dfftw_execute_dft_c2r(back_plan, dfdy_out(:,:,k), dfdy_out(:,:,k))
end do

end subroutine filt_da_cpu

subroutine filt_da_gpu_planmany(field_in, f_out, dfdx_out, dfdy_out, kx_in,    &
    ky_in)
implicit none
real(dp), intent(in) :: field_in(ld, ny, nbatch)
real(dp), intent(out) :: f_out(ld, ny, nbatch)
real(dp), intent(out) :: dfdx_out(ld, ny, nbatch)
real(dp), intent(out) :: dfdy_out(ld, ny, nbatch)
real(dp), intent(in) :: kx_in(lh, ny)
real(dp), intent(in) :: ky_in(lh, ny)

integer :: fw_plan, bk_plan, istat
integer :: n_s(2), inem_s(2), onem_s(2)
integer :: i, j, k
real(dp) :: const, fr, fi
real(dp), device, allocatable :: f_d(:,:,:)
real(dp), device, allocatable :: dfdx_d(:,:,:)
real(dp), device, allocatable :: dfdy_d(:,:,:)
real(dp), device, allocatable :: kx_d(:,:)
real(dp), device, allocatable :: ky_d(:,:)
complex(dp), device, allocatable :: fh_d(:,:,:)
complex(dp), device, allocatable :: dfdxh_d(:,:,:)
complex(dp), device, allocatable :: dfdyh_d(:,:,:)

const = 1.0_dp / real(nx * ny, dp)
allocate(f_d(ld, ny, nbatch), dfdx_d(ld, ny, nbatch),                         &
    dfdy_d(ld, ny, nbatch))
allocate(kx_d(lh, ny), ky_d(lh, ny))
allocate(fh_d(lh, ny, nbatch), dfdxh_d(lh, ny, nbatch),                       &
    dfdyh_d(lh, ny, nbatch))

kx_d = kx_in
ky_d = ky_in
f_d = const * field_in

n_s = (/ ny, nx /)
inem_s = (/ ny, ld /)
onem_s = (/ ny, lh /)

istat = cufftPlanMany(fw_plan, 2, n_s, inem_s, 1, ld*ny, onem_s, 1,           &
    lh*ny, CUFFT_D2Z, nbatch)
call require_cufft_success('forward plan', istat)

istat = cufftPlanMany(bk_plan, 2, n_s, onem_s, 1, lh*ny, inem_s, 1,           &
    ld*ny, CUFFT_Z2D, nbatch)
call require_cufft_success('inverse plan', istat)

istat = cufftExecD2Z(fw_plan, f_d, fh_d)
call require_cufft_success('forward exec', istat)
istat = cudaDeviceSynchronize()
call require_cuda_success('forward sync', istat)

!$cuf kernel do(3) <<<*,*>>>
do k = 1, nbatch
do j = 1, ny
do i = 1, lh
    if (i == lh .or. j == ny/2 + 1) then
        fh_d(i,j,k) = cmplx(0.0_dp, 0.0_dp, kind=dp)
        dfdxh_d(i,j,k) = cmplx(0.0_dp, 0.0_dp, kind=dp)
        dfdyh_d(i,j,k) = cmplx(0.0_dp, 0.0_dp, kind=dp)
    else
        fr = real(fh_d(i,j,k), kind=dp)
        fi = aimag(fh_d(i,j,k))
        dfdxh_d(i,j,k) = cmplx(-fi * kx_d(i,j), fr * kx_d(i,j), kind=dp)
        dfdyh_d(i,j,k) = cmplx(-fi * ky_d(i,j), fr * ky_d(i,j), kind=dp)
    end if
end do
end do
end do

istat = cufftExecZ2D(bk_plan, fh_d, f_d)
call require_cufft_success('inverse filtered exec', istat)
istat = cufftExecZ2D(bk_plan, dfdxh_d, dfdx_d)
call require_cufft_success('inverse dfdx exec', istat)
istat = cufftExecZ2D(bk_plan, dfdyh_d, dfdy_d)
call require_cufft_success('inverse dfdy exec', istat)
istat = cudaDeviceSynchronize()
call require_cuda_success('inverse sync', istat)

f_out = f_d
dfdx_out = dfdx_d
dfdy_out = dfdy_d

istat = cufftDestroy(fw_plan)
call require_cufft_success('destroy forward plan', istat)
istat = cufftDestroy(bk_plan)
call require_cufft_success('destroy inverse plan', istat)

end subroutine filt_da_gpu_planmany

subroutine check_array(name, expected, actual, worst)
implicit none
character(len=*), intent(in) :: name
real(dp), intent(in) :: expected(ld, ny, nbatch)
real(dp), intent(in) :: actual(ld, ny, nbatch)
real(dp), intent(inout) :: worst
real(dp) :: err, ref, rel

err = maxval(abs(expected - actual))
ref = maxval(abs(expected))
rel = err / max(ref, 1.0_dp)
worst = max(worst, err)

write(*,'(a,2x,a,es12.4,2x,a,es12.4)') trim(name), 'max_abs=', err,          &
    'rel=', rel

end subroutine check_array

subroutine require_cufft_success(where, istat)
implicit none
character(len=*), intent(in) :: where
integer, intent(in) :: istat

if (istat /= CUFFT_SUCCESS) then
    write(*,'(a,a,i0)') 'cuFFT failure at ', trim(where), istat
    stop 2
end if

end subroutine require_cufft_success

subroutine require_cuda_success(where, istat)
implicit none
character(len=*), intent(in) :: where
integer, intent(in) :: istat

if (istat /= 0) then
    write(*,'(a,a,i0)') 'CUDA failure at ', trim(where), istat
    stop 3
end if

end subroutine require_cuda_success

end program validate_filt_da_cufft
