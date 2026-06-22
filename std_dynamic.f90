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
subroutine std_dynamic(Cs_1D)
!*******************************************************************************
!
! Subroutine uses the standard dynamic model to calculate the Smagorinsky
! coefficient Cs_1D and |S|. This is done layer-by-layer to save memory.
!
use types, only : rprec
use param, only : ld, ny, nz, coord
use test_filtermodule
use sim_param, only : u, v, w
use sgs_param, only : ee_now, S11, S12, S13, S22, S23, S33, delta, S,          &
    u_bar, v_bar, w_bar, L11, L12, L13, L22, L23, L33,                         &
    M11, M12, M13, M22, M23, M33,                                              &
    S_bar, S11_bar, S12_bar,S13_bar, S22_bar, S23_bar, S33_bar,                &
    S_S11_bar, S_S12_bar, S_S13_bar, S_S22_bar, S_S23_bar, S_S33_bar
#ifdef ENABLE_CUDA
use cudafor
#endif
implicit none

real(rprec), dimension(nz), intent(out) :: Cs_1D
real(rprec), dimension(ld,ny) :: LM, MM
real(rprec) :: const, lm_sum, mm_sum, cs_val
integer :: jx, jy, jz
#ifdef ENABLE_CUDA
logical :: std_dynamic_cuda_enabled
logical :: use_std_cuda

use_std_cuda = std_dynamic_cuda_enabled()
#endif

do jz = 1, nz
    ! using L_ij as temp storage here
    ! watch the 0.25's:  recall w = c*z^2 close to wall, so get 0.25
#ifdef ENABLE_CUDA
    if (use_std_cuda) then
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            if ((coord == 0) .and. (jz == 1)) then
                L11(jx,jy) = u(jx,jy,1)*u(jx,jy,1)
                L12(jx,jy) = u(jx,jy,1)*v(jx,jy,1)
                L13(jx,jy) = u(jx,jy,1)*0.25_rprec*w(jx,jy,2)
                L22(jx,jy) = v(jx,jy,1)*v(jx,jy,1)
                L23(jx,jy) = v(jx,jy,jz)*0.25_rprec*w(jx,jy,2)
                L33(jx,jy) = (0.25_rprec*w(jx,jy,2))**2
                u_bar(jx,jy) = u(jx,jy,1)
                v_bar(jx,jy) = v(jx,jy,1)
                w_bar(jx,jy) = 0.25_rprec*w(jx,jy,2)
            else
                u_bar(jx,jy) = 0.5_rprec*(u(jx,jy,jz) + u(jx,jy,jz-1))
                v_bar(jx,jy) = 0.5_rprec*(v(jx,jy,jz) + v(jx,jy,jz-1))
                w_bar(jx,jy) = w(jx,jy,jz)
                L11(jx,jy) = u_bar(jx,jy)*u_bar(jx,jy)
                L12(jx,jy) = u_bar(jx,jy)*v_bar(jx,jy)
                L13(jx,jy) = u_bar(jx,jy)*w_bar(jx,jy)
                L22(jx,jy) = v_bar(jx,jy)*v_bar(jx,jy)
                L23(jx,jy) = v_bar(jx,jy)*w_bar(jx,jy)
                L33(jx,jy) = w_bar(jx,jy)*w_bar(jx,jy)
            end if
        end do
        end do
        call std_dynamic_cuda_sync('velocity and L setup')
    else
#endif
    if ( (coord == 0) .and. (jz == 1) ) then
        ! put on uvp-nodes
        L11(:,:) = u(:,:,1)*u(:,:,1)                ! uv-node
        L12(:,:) = u(:,:,1)*v(:,:,1)                ! uv-node
        L13(:,:) = u(:,:,1)*0.25_rprec*w(:,:,2)     ! parabolic interp.
        L22(:,:) = v(:,:,1)*v(:,:,1)                ! uv-node
        L23(:,:) = v(:,:,jz)*0.25_rprec*w(:,:,2)    ! uv-node
        L33(:,:) = (0.25_rprec*w(:,:,2))**2         ! uv-node
        u_bar(:,:) = u(:,:,1)
        v_bar(:,:) = v(:,:,1)
        w_bar(:,:) = 0.25_rprec*w(:,:,2)
    else
        ! w-nodes
        L11(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))*                        &
            0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))
        L12(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))*                        &
            0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
        L13(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))*w(:,:,jz)
        L22(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))*                        &
            0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
        L23(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))*w(:,:,jz)
        L33(:,:) = w(:,:,jz)*w(:,:,jz)
        u_bar(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))
        v_bar(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
        w_bar(:,:) = w(:,:,jz)
   end if
#ifdef ENABLE_CUDA
    end if
#endif

    ! in-place filtering
    call test_filter_3 ( u_bar, v_bar, w_bar )

#ifdef ENABLE_CUDA
    if (use_std_cuda) then
        call test_filter_6 ( L11, L12, L13, L22, L23, L33 )
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            L11(jx,jy) = L11(jx,jy) - u_bar(jx,jy)*u_bar(jx,jy)
            L12(jx,jy) = L12(jx,jy) - u_bar(jx,jy)*v_bar(jx,jy)
            L13(jx,jy) = L13(jx,jy) - u_bar(jx,jy)*w_bar(jx,jy)
            L22(jx,jy) = L22(jx,jy) - v_bar(jx,jy)*v_bar(jx,jy)
            L23(jx,jy) = L23(jx,jy) - v_bar(jx,jy)*w_bar(jx,jy)
            L33(jx,jy) = L33(jx,jy) - w_bar(jx,jy)*w_bar(jx,jy)
        end do
        end do
        call std_dynamic_cuda_sync('filtered Lij correction')
    else
#endif
    call test_filter ( L11 )
    L11 = L11 - u_bar*u_bar
    call test_filter ( L12 )
    L12 = L12 - u_bar*v_bar
    call test_filter ( L13 )
    L13 = L13 - u_bar*w_bar
    call test_filter ( L22 )
    L22 = L22 - v_bar*v_bar
    call test_filter ( L23 )
    L23 = L23 - v_bar*w_bar
    call test_filter ( L33 )
    L33 = L33 - w_bar*w_bar
#ifdef ENABLE_CUDA
    end if
#endif

    ! calculate |S|
#ifdef ENABLE_CUDA
    if (use_std_cuda) then
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            S(jx,jy) = sqrt(2._rprec*(S11(jx,jy,jz)**2 + S22(jx,jy,jz)**2 +      &
                S33(jx,jy,jz)**2 + 2._rprec*(S12(jx,jy,jz)**2 +                 &
                S13(jx,jy,jz)**2 + S23(jx,jy,jz)**2)))
            S11_bar(jx,jy) = S11(jx,jy,jz)
            S12_bar(jx,jy) = S12(jx,jy,jz)
            S13_bar(jx,jy) = S13(jx,jy,jz)
            S22_bar(jx,jy) = S22(jx,jy,jz)
            S23_bar(jx,jy) = S23(jx,jy,jz)
            S33_bar(jx,jy) = S33(jx,jy,jz)
            S_S11_bar(jx,jy) = S(jx,jy)*S11(jx,jy,jz)
            S_S12_bar(jx,jy) = S(jx,jy)*S12(jx,jy,jz)
            S_S13_bar(jx,jy) = S(jx,jy)*S13(jx,jy,jz)
            S_S22_bar(jx,jy) = S(jx,jy)*S22(jx,jy,jz)
            S_S23_bar(jx,jy) = S(jx,jy)*S23(jx,jy,jz)
            S_S33_bar(jx,jy) = S(jx,jy)*S33(jx,jy,jz)
        end do
        end do
        call std_dynamic_cuda_sync('S setup')
        call test_filter_6 ( S11_bar, S12_bar, S13_bar, S22_bar, S23_bar,      &
            S33_bar )
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            S_bar(jx,jy) = sqrt(2._rprec*(S11_bar(jx,jy)**2 +                 &
                S22_bar(jx,jy)**2 + S33_bar(jx,jy)**2 + 2._rprec*(            &
                S12_bar(jx,jy)**2 + S13_bar(jx,jy)**2 + S23_bar(jx,jy)**2)))
        end do
        end do
        call std_dynamic_cuda_sync('filtered S magnitude')
    else
#endif
    S(:,:) = sqrt(2._rprec*(S11(:,:,jz)**2 + S22(:,:,jz)**2 +                  &
        S33(:,:,jz)**2 + 2._rprec*(S12(:,:,jz)**2 +                            &
        S13(:,:,jz)**2 + S23(:,:,jz)**2)))

    ! S_ij already on w-nodes
    S11_bar(:,:) = S11(:,:,jz)
    S12_bar(:,:) = S12(:,:,jz)
    S13_bar(:,:) = S13(:,:,jz)
    S22_bar(:,:) = S22(:,:,jz)
    S23_bar(:,:) = S23(:,:,jz)
    S33_bar(:,:) = S33(:,:,jz)

    call test_filter ( S11_bar )
    call test_filter ( S12_bar )
    call test_filter ( S13_bar )
    call test_filter ( S22_bar )
    call test_filter ( S23_bar )
    call test_filter ( S33_bar )

    S_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 +              &
        2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))

    S_S11_bar(:,:) = S(:,:)*S11(:,:,jz)
    S_S12_bar(:,:) = S(:,:)*S12(:,:,jz)
    S_S13_bar(:,:) = S(:,:)*S13(:,:,jz)
    S_S22_bar(:,:) = S(:,:)*S22(:,:,jz)
    S_S23_bar(:,:) = S(:,:)*S23(:,:,jz)
    S_S33_bar(:,:) = S(:,:)*S33(:,:,jz)
#ifdef ENABLE_CUDA
    end if
#endif

    call test_filter_6 ( S_S11_bar, S_S12_bar, S_S13_bar, S_S22_bar,           &
        S_S23_bar, S_S33_bar )

    ! now put beta back into M_ij
    const = 2._rprec*delta**2
#ifdef ENABLE_CUDA
    if (use_std_cuda) then
        lm_sum = 0._rprec
        mm_sum = 0._rprec
!$cuf kernel do(2) <<<*,*>>> reduction(+:lm_sum,mm_sum)
        do jy = 1, ny
        do jx = 1, ld
            M11(jx,jy) = const*(S_S11_bar(jx,jy) -                              &
                4._rprec*S_bar(jx,jy)*S11_bar(jx,jy))
            M12(jx,jy) = const*(S_S12_bar(jx,jy) -                              &
                4._rprec*S_bar(jx,jy)*S12_bar(jx,jy))
            M13(jx,jy) = const*(S_S13_bar(jx,jy) -                              &
                4._rprec*S_bar(jx,jy)*S13_bar(jx,jy))
            M22(jx,jy) = const*(S_S22_bar(jx,jy) -                              &
                4._rprec*S_bar(jx,jy)*S22_bar(jx,jy))
            M23(jx,jy) = const*(S_S23_bar(jx,jy) -                              &
                4._rprec*S_bar(jx,jy)*S23_bar(jx,jy))
            M33(jx,jy) = const*(S_S33_bar(jx,jy) -                              &
                4._rprec*S_bar(jx,jy)*S33_bar(jx,jy))
            lm_sum = lm_sum + L11(jx,jy)*M11(jx,jy) + L22(jx,jy)*M22(jx,jy) +    &
                L33(jx,jy)*M33(jx,jy) + 2._rprec*(L12(jx,jy)*M12(jx,jy) +        &
                L13(jx,jy)*M13(jx,jy) + L23(jx,jy)*M23(jx,jy))
            mm_sum = mm_sum + M11(jx,jy)**2 + M22(jx,jy)**2 + M33(jx,jy)**2 +    &
                2._rprec*(M12(jx,jy)**2 + M13(jx,jy)**2 + M23(jx,jy)**2)
        end do
        end do
        call std_dynamic_cuda_sync('Mij reductions')
        cs_val = max(0._rprec, lm_sum/mm_sum)
        Cs_1D(jz) = cs_val
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            ee_now(jx,jy,jz) = L11(jx,jy)**2 + L22(jx,jy)**2 +                  &
                L33(jx,jy)**2 + 2._rprec*(L12(jx,jy)**2 + L13(jx,jy)**2 +        &
                L23(jx,jy)**2) - 2._rprec*(L11(jx,jy)*M11(jx,jy) +              &
                L22(jx,jy)*M22(jx,jy) + L33(jx,jy)*M33(jx,jy) +                 &
                2._rprec*(L12(jx,jy)*M12(jx,jy) + L13(jx,jy)*M13(jx,jy) +        &
                L23(jx,jy)*M23(jx,jy)))*cs_val + (M11(jx,jy)**2 +               &
                M22(jx,jy)**2 + M33(jx,jy)**2 + 2._rprec*(M12(jx,jy)**2 +       &
                M13(jx,jy)**2 + M23(jx,jy)**2))*cs_val**2
        end do
        end do
        call std_dynamic_cuda_sync('ee_now')
    else
#endif
    M11 = const*(S_S11_bar - 4._rprec*S_bar*S11_bar)
    M12 = const*(S_S12_bar - 4._rprec*S_bar*S12_bar)
    M13 = const*(S_S13_bar - 4._rprec*S_bar*S13_bar)
    M22 = const*(S_S22_bar - 4._rprec*S_bar*S22_bar)
    M23 = const*(S_S23_bar - 4._rprec*S_bar*S23_bar)
    M33 = const*(S_S33_bar - 4._rprec*S_bar*S33_bar)

    Cs_1D(jz) =                                                                &
        sum(L11*M11+L22*M22+L33*M33+2._rprec*(L12*M12+L13*M13+L23*M23))/       &
        sum(M11**2 + M22**2 + M33**2 + 2._rprec*(M12**2 + M13**2 + M23**2))
    Cs_1D(jz) = max(0._rprec, Cs_1D(jz))

    ! Calculate ee_now (the current value of eij*eij)
    LM = L11*M11+L22*M22+L33*M33+2._rprec*(L12*M12+L13*M13+L23*M23)
    MM = M11**2+M22**2+M33**2+2._rprec*(M12**2+M13**2+M23**2)
    ee_now(:,:,jz) = L11**2+L22**2+L33**2+2._rprec*(L12**2+L13**2+L23**2)      &
                     -2._rprec*LM*Cs_1D(jz) + MM*Cs_1D(jz)**2
#ifdef ENABLE_CUDA
    end if
#endif
end do

end subroutine std_dynamic

#ifdef ENABLE_CUDA
!*******************************************************************************
logical function std_dynamic_cuda_enabled()
!*******************************************************************************
implicit none

std_dynamic_cuda_enabled = .true.

end function std_dynamic_cuda_enabled

!*******************************************************************************
subroutine std_dynamic_cuda_sync(where)
!*******************************************************************************
use cudafor
implicit none

character(len=*), intent(in) :: where
integer :: istat

istat = cudaDeviceSynchronize()
if (istat /= 0) then
    print *, 'std_dynamic CUDA sync failure at ', trim(where), ': ', istat
    stop
end if
istat = cudaGetLastError()
if (istat /= 0) then
    print *, 'std_dynamic CUDA kernel failure at ', trim(where), ': ', istat
    stop
end if

end subroutine std_dynamic_cuda_sync
#endif
