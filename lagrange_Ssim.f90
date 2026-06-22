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
subroutine lagrange_Ssim()
!*******************************************************************************
!
! This subroutine dynamically calculates Cs_opt2. See  Meneveau, Lund, Cabot,
! JFM, 319: 353-385 (1996) DOI: 10.1017/S0022112096007379
!
use types, only : rprec
use param
use sim_param, only : u, v, w
use sgs_param, only : F_LM, F_MM, Beta, Cs_opt2, opftime, lagran_dt
use sgs_param, only : S11, S12, S13, S22, S23, S33, delta, S, ee_now, Tn_all,  &
    u_bar, v_bar, w_bar,                                                       &
    L11, L12, L13, L22, L23, L33, M11, M12, M13, M22, M23, M33,                &
    S_bar, S11_bar, S12_bar, S13_bar, S22_bar, S23_bar, S33_bar,               &
    S_S11_bar, S_S12_bar, S_S13_bar, S_S22_bar, S_S23_bar, S_S33_bar
use test_filtermodule
use messages
use string_util, only : string_concat
#ifdef PPLVLSET
use level_set, only : level_set_lag_dyn, level_set_Cs_lag_dyn
#endif
#ifdef PPDYN_TN
use sgs_param, only:F_ee2,F_deedt2,ee_past
#endif
#ifdef PPMPI
use mpi_defs, only:mpi_sync_real_array,MPI_SYNC_DOWNUP
#endif
#ifdef ENABLE_CUDA
use cudafor
#endif

implicit none

character(*), parameter :: sub_name = 'lagrange_Ssim'

real(rprec), parameter :: eps = 1.e-32_rprec
#ifdef ENABLE_CUDA
real(rprec), managed, allocatable, save, dimension(:,:) :: fourbeta
real(rprec), managed, allocatable, save, dimension(:,:) :: LM, MM, Tn, epsi, dumfac
#else
real(rprec), dimension(ld,ny) :: fourbeta
real(rprec), dimension(ld,ny) :: LM, MM, Tn, epsi, dumfac
#endif
real(rprec) :: const
real(rprec) :: opftdelta
integer :: istart, iend, iend_mod, jx, jy, jz, ii, i
logical, save :: F_LM_MM_init = .false.
#ifdef ENABLE_CUDA
logical :: lagrange_ssim_cuda_enabled
logical :: use_lag_ssim_cuda

if (.not. allocated(fourbeta)) then
    allocate(fourbeta(ld,ny), LM(ld,ny), MM(ld,ny), Tn(ld,ny), epsi(ld,ny),      &
        dumfac(ld,ny))
end if
use_lag_ssim_cuda = lagrange_ssim_cuda_enabled()
#endif

! Set coefficients
opftdelta = opftime*delta
const = 2._rprec*delta*delta

#ifdef PPLVLSET
call level_set_lag_dyn (S11, S12, S13, S22, S23, S33)
#else
#ifdef ENABLE_CUDA
if (use_lag_ssim_cuda) then
!$cuf kernel do(3) <<<*,*>>>
    do jz = 1, nz
    do jy = 1, ny
    do jx = 1, ld
        Beta(jx,jy,jz) = 1._rprec
    end do
    end do
    end do
    call lagrange_ssim_cuda_sync('initialize Beta')
else
#endif
Beta = 1._rprec
#ifdef ENABLE_CUDA
end if
#endif
#endif

! "Rearrange" F_LM, F_MM, F_ee2, F_deedt2 (running averages) so that
!   their new positions (i,j,k) correspond to the current (i,j,k) particle
call interpolag_Ssim()

! For each horizontal level, calculate Lij(:,:) and Mij(:,:).  Then update
!   the running averages, F_LM(:,:,jz) and F_MM(:,:,jz), which are used to
!   calculate Cs_opt2(:,:,jz).
do jz = 1,nz

    ! Calculate Lij
    ! Interp u,v,w onto w-nodes and store result as u_bar,v_bar,w_bar
    ! (except for very first level which should be on uvp-nodes)
#ifdef ENABLE_CUDA
    if (use_lag_ssim_cuda) then
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            if ((coord == 0) .and. (jz == 1)) then
                u_bar(jx,jy) = u(jx,jy,1)
                v_bar(jx,jy) = v(jx,jy,1)
                w_bar(jx,jy) = 0.25_rprec*w(jx,jy,2)
            else
                u_bar(jx,jy) = 0.5_rprec*(u(jx,jy,jz) + u(jx,jy,jz-1))
                v_bar(jx,jy) = 0.5_rprec*(v(jx,jy,jz) + v(jx,jy,jz-1))
                w_bar(jx,jy) = w(jx,jy,jz)
            end if
            L11(jx,jy) = u_bar(jx,jy)*u_bar(jx,jy)
            L12(jx,jy) = u_bar(jx,jy)*v_bar(jx,jy)
            L13(jx,jy) = u_bar(jx,jy)*w_bar(jx,jy)
            L23(jx,jy) = v_bar(jx,jy)*w_bar(jx,jy)
            L22(jx,jy) = v_bar(jx,jy)*v_bar(jx,jy)
            L33(jx,jy) = w_bar(jx,jy)*w_bar(jx,jy)
        end do
        end do
        call lagrange_ssim_cuda_sync('velocity and L setup')
    else
#endif
    if ( ( coord == 0 ) .and. (jz == 1) ) then  ! uvp-nodes
        u_bar(:,:) = u(:,:,1)
        v_bar(:,:) = v(:,:,1)
        w_bar(:,:) = .25_rprec*w(:,:,2)
    else  ! w-nodes
        u_bar(:,:) = .5_rprec*(u(:,:,jz) + u(:,:,jz-1))
        v_bar(:,:) = .5_rprec*(v(:,:,jz) + v(:,:,jz-1))
        w_bar(:,:) = w(:,:,jz)
    end if

    ! First term before filtering (not the final value)
    L11 = u_bar*u_bar
    L12 = u_bar*v_bar
    L13 = u_bar*w_bar
    L23 = v_bar*w_bar
    L22 = v_bar*v_bar
    L33 = w_bar*w_bar
#ifdef ENABLE_CUDA
    end if
#endif

    ! Filter first term and add the second term to get the final value
    call test_filter_3 ( u_bar, v_bar, w_bar )   ! in-place filtering
#ifdef ENABLE_CUDA
    if (use_lag_ssim_cuda) then
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
        call lagrange_ssim_cuda_sync('filtered Lij correction')
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

    ! Calculate |S|
#ifdef ENABLE_CUDA
    if (use_lag_ssim_cuda) then
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            S(jx,jy) = sqrt(2._rprec*(S11(jx,jy,jz)**2 + S22(jx,jy,jz)**2 +    &
                S33(jx,jy,jz)**2 + 2._rprec*(S12(jx,jy,jz)**2 +               &
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
        call lagrange_ssim_cuda_sync('S setup')
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
        call lagrange_ssim_cuda_sync('filtered S magnitude')
    else
#endif
    S(:,:) = sqrt(2._rprec*(S11(:,:,jz)**2+S22(:,:,jz)**2+S33(:,:,jz)**2+&
        2._rprec*(S12(:,:,jz)**2+S13(:,:,jz)**2+S23(:,:,jz)**2)))

    ! Select Sij for this level for test-filtering, saving results as Sij_bar
    !   note: Sij is already on w-nodes
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

    ! Calculate |S_bar| (the test-filtered Sij)
    S_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 +&
        2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))

    ! Calculate |S|Sij then test-filter this quantity
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

    ! Calculate Mij
#ifdef ENABLE_CUDA
    if (use_lag_ssim_cuda) then
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            fourbeta(jx,jy) = 4._rprec*Beta(jx,jy,jz)
            M11(jx,jy) = const*(S_S11_bar(jx,jy) -                             &
                fourbeta(jx,jy)*S_bar(jx,jy)*S11_bar(jx,jy))
            M12(jx,jy) = const*(S_S12_bar(jx,jy) -                             &
                fourbeta(jx,jy)*S_bar(jx,jy)*S12_bar(jx,jy))
            M13(jx,jy) = const*(S_S13_bar(jx,jy) -                             &
                fourbeta(jx,jy)*S_bar(jx,jy)*S13_bar(jx,jy))
            M22(jx,jy) = const*(S_S22_bar(jx,jy) -                             &
                fourbeta(jx,jy)*S_bar(jx,jy)*S22_bar(jx,jy))
            M23(jx,jy) = const*(S_S23_bar(jx,jy) -                             &
                fourbeta(jx,jy)*S_bar(jx,jy)*S23_bar(jx,jy))
            M33(jx,jy) = const*(S_S33_bar(jx,jy) -                             &
                fourbeta(jx,jy)*S_bar(jx,jy)*S33_bar(jx,jy))
            LM(jx,jy) = L11(jx,jy)*M11(jx,jy) + L22(jx,jy)*M22(jx,jy) +        &
                L33(jx,jy)*M33(jx,jy) + 2._rprec*(L12(jx,jy)*M12(jx,jy) +      &
                L13(jx,jy)*M13(jx,jy) + L23(jx,jy)*M23(jx,jy))
            MM(jx,jy) = M11(jx,jy)**2 + M22(jx,jy)**2 + M33(jx,jy)**2 +        &
                2._rprec*(M12(jx,jy)**2 + M13(jx,jy)**2 + M23(jx,jy)**2)
            ee_now(jx,jy,jz) = L11(jx,jy)**2 + L22(jx,jy)**2 +                &
                L33(jx,jy)**2 + 2._rprec*(L12(jx,jy)**2 + L13(jx,jy)**2 +      &
                L23(jx,jy)**2) - 2._rprec*LM(jx,jy)*Cs_opt2(jx,jy,jz) +        &
                MM(jx,jy)*Cs_opt2(jx,jy,jz)**2
        end do
        end do
        call lagrange_ssim_cuda_sync('M contractions')
    else
#endif
    fourbeta=4._rprec*Beta(:,:,jz)

    M11 = const*(S_S11_bar - fourbeta*S_bar*S11_bar)
    M12 = const*(S_S12_bar - fourbeta*S_bar*S12_bar)
    M13 = const*(S_S13_bar - fourbeta*S_bar*S13_bar)
    M22 = const*(S_S22_bar - fourbeta*S_bar*S22_bar)
    M23 = const*(S_S23_bar - fourbeta*S_bar*S23_bar)
    M33 = const*(S_S33_bar - fourbeta*S_bar*S33_bar)

    ! Calculate LijMij and MijMij for each point in the plane
    LM = L11*M11+L22*M22+L33*M33+2._rprec*(L12*M12+L13*M13+L23*M23)
    MM = M11**2+M22**2+M33**2+2._rprec*(M12**2+M13**2+M23**2)

    ! Calculate ee_now (the current value of eij*eij)
    ee_now(:,:,jz) = L11**2+L22**2+L33**2+2._rprec*(L12**2+L13**2+L23**2)      &
        -2._rprec*LM*Cs_opt2(:,:,jz) + MM*Cs_opt2(:,:,jz)**2
#ifdef ENABLE_CUDA
    end if
#endif

    ! Using local time counter to reinitialize SGS quantities when restarting
    if (inilag) then
        if ((.not. F_LM_MM_init) .and. (jt == cs_count .or. jt == DYN_init)) then
            print *,'F_MM and F_LM initialized'
#ifdef ENABLE_CUDA
            if (use_lag_ssim_cuda) then
!$cuf kernel do(2) <<<*,*>>>
                do jy = 1, ny
                do jx = 1, ld
                    if (jx >= ld-1) then
                        F_MM(jx,jy,jz) = 1._rprec
                        F_LM(jx,jy,jz) = 1._rprec
                    else
                        F_MM(jx,jy,jz) = MM(jx,jy)
                        F_LM(jx,jy,jz) = 0.025_rprec*MM(jx,jy)
                    end if
                end do
                end do
                call lagrange_ssim_cuda_sync('initialize F_LM/F_MM')
            else
#endif
            F_MM (:,:,jz) = MM
            F_LM (:,:,jz) = 0.025_rprec*MM
            F_MM(ld-1:ld,:,jz)=1._rprec
            F_LM(ld-1:ld,:,jz)=1._rprec
#ifdef ENABLE_CUDA
            end if
#endif

            if (jz == 1) then
                if (coord == 0) then
                    write(*, *) 'LM(1, 1)=', LM(1, 1)
                    write(*, *) 'MM(1, 1)=', MM(1, 1)
                    write(*, *) 'M11(1, 1)=', M11(1, 1)
                    write(*, *) 'S_S11_bar(1, 1)=', S_S11_bar(1, 1)
                    write(*, *) 'S11(1, 1, 1)=', S11(1, 1, 1)
                    write(*, *) 'S(1, 1)=', S(1, 1)
                    write(*, *) 'S11_bar(1, 1)=', S11_bar(1, 1)
                    write(*, *) 'S_bar(1, 1)=', S_bar(1, 1)
                endif
            endif

            if (jz == nz) F_LM_MM_init = .true.
        endif
    endif

    ! Inflow
    if (inflow_type > 0) then
        iend = floor (fringe_region_end * nx + 1._rprec)
        istart = floor ((fringe_region_end - fringe_region_len) * nx + 1._rprec)

#ifdef ENABLE_CUDA
        if (use_lag_ssim_cuda) then
            if (istart > iend) then
                write (*, *) 'lagrange_Ssim: istart > iend'
                stop
            endif
            iend_mod = modulo(iend - 1, nx) + 1
!$cuf kernel do(2) <<<*,*>>>
            do jy = 1, ny
            do jx = 1, ld
                Tn(jx,jy) = max(MM(jx,jy), 0.1_rprec*const*S(jx,jy)**2)
                MM(jx,jy) = Tn(jx,jy)
                if (iend <= nx) then
                    if ((jx >= istart) .and. (jx <= iend)) then
                        LM(jx,jy) = 0._rprec
                        F_LM(jx,jy,jz) = 0._rprec
                    end if
                else
                    if ((jx >= istart) .or. (jx <= iend_mod)) then
                        LM(jx,jy) = 0._rprec
                        F_LM(jx,jy,jz) = 0._rprec
                    end if
                end if
            end do
            end do
            call lagrange_ssim_cuda_sync('inflow SGS damping')
        else
#endif
        Tn = merge(.1_rprec*const*S**2,MM,MM.le..1_rprec*const*S**2)
        MM = Tn

        if (istart > iend) then
            write (*, *) 'lagrange_Ssim: istart > iend'
            stop
        endif
        do i = istart, iend
            ii = modulo (i - 1, nx) + 1
            LM(ii, :) = 0._rprec
            F_LM(ii, :, jz) = 0._rprec
        enddo
#ifdef ENABLE_CUDA
        end if
#endif
    endif

    ! Update running averages (F_LM, F_MM, F_ee2, F_deedt2)
    ! Determine averaging timescale
#ifdef ENABLE_CUDA
#ifndef PPDYN_TN
    if (use_lag_ssim_cuda) then
!$cuf kernel do(2) <<<*,*>>>
        do jy = 1, ny
        do jx = 1, ld
            Tn(jx,jy) = max(F_LM(jx,jy,jz)*F_MM(jx,jy,jz), eps)
            Tn(jx,jy) = opftdelta/sqrt(sqrt(sqrt(Tn(jx,jy))))
            dumfac(jx,jy) = lagran_dt/Tn(jx,jy)
            epsi(jx,jy) = dumfac(jx,jy)/(1._rprec + dumfac(jx,jy))
            F_LM(jx,jy,jz) = epsi(jx,jy)*LM(jx,jy) +                            &
                (1._rprec-epsi(jx,jy))*F_LM(jx,jy,jz)
            F_MM(jx,jy,jz) = epsi(jx,jy)*MM(jx,jy) +                            &
                (1._rprec-epsi(jx,jy))*F_MM(jx,jy,jz)
            F_LM(jx,jy,jz) = max(eps, F_LM(jx,jy,jz))
            Cs_opt2(jx,jy,jz) = F_LM(jx,jy,jz)/(F_MM(jx,jy,jz) + eps)
            if (jx >= ld-1) Cs_opt2(jx,jy,jz) = 0._rprec
            Cs_opt2(jx,jy,jz) = max(eps, Cs_opt2(jx,jy,jz))
            Tn_all(jx,jy,jz) = Tn(jx,jy)
        end do
        end do
        call lagrange_ssim_cuda_sync('update F_LM/F_MM and Cs')
    else
#endif
#endif
#ifdef PPDYN_TN
    ! based on Taylor timescale
    Tn = 4._rprec*pi*sqrt(F_ee2(:,:,jz)/F_deedt2(:,:,jz))
#else
    ! based on Meneveau, Cabot, Lund paper (JFM 1996)
    Tn = max (F_LM(:,:,jz) * F_MM(:,:,jz), eps)
    Tn = opftdelta/sqrt(sqrt(sqrt(Tn)))
#endif

    ! Calculate new running average = old*(1-epsi) + instantaneous*epsi
    dumfac = lagran_dt/Tn
    epsi = dumfac / (1._rprec + dumfac)

    F_LM(:,:,jz)=(epsi*LM + (1.0_rprec-epsi)*F_LM(:,:,jz))
    F_MM(:,:,jz)=(epsi*MM + (1.0_rprec-epsi)*F_MM(:,:,jz))
    ! clipping to avoid instability
    F_LM(:,:,jz)= max(eps, F_LM(:,:,jz))

#ifdef PPDYN_TN
    ! note: the instantaneous value of the derivative is a Lagrangian average
    F_ee2(:,:,jz) = epsi*ee_now(:,:,jz)**2 + (1._rprec-epsi)*F_ee2(:,:,jz)
    F_deedt2(:,:,jz) = epsi*( ((ee_now(:,:,jz)-ee_past(:,:,jz))/lagran_dt)**2 )&
        + (1._rprec-epsi)*F_deedt2(:,:,jz)
    ee_past(:,:,jz) = ee_now(:,:,jz)
#endif

    ! Calculate Cs_opt2 (use only one of the methods below)
    ! Standard method - LASS
    ! Added +eps to avoid division by zero
    Cs_opt2(:,:,jz) = F_LM(:,:,jz) / (F_MM(:,:,jz) + eps)
    Cs_opt2(ld-1:ld,:,jz) = 0._rprec

    ! 9-point average
    !do i=1,nx
    !do j=1,ny
    !    ilo=i-1; ihi=i+1; jlo=j-1;  jhi=j+1
    !    if (ilo.eq.0) ilo=nx
    !    if (jlo.eq.0) jlo=ny
    !    if (ihi.eq.nx+1) ihi=1
    !    if (jhi.eq.ny+1) jhi=1
    !    Cs_opt2(i,j,jz) = (LM(i,j)+LM(ilo,j)+LM(ihi,j)+LM(ilo,jlo)+&
    !                       LM(ihi,jlo)+LM(i,jlo)+LM(ilo,jhi)+LM(i,jhi)+LM(ihi,jhi))/ &
    !                      (MM(i,j)+MM(ilo,j)+MM(ihi,j)+MM(ilo,jlo)+&
    !                       MM(ihi,jlo)+MM(i,jlo)+MM(ilo,jhi)+MM(i,jhi)+MM(ihi,jhi))
    !enddo
    !enddo

    ! Directly
    !Cs_opt2(:,:,jz) = LM(:,:)/MM(:,:)

    ! Clip Cs if necessary
    Cs_opt2(:,:,jz)= max (eps, Cs_opt2(:,:,jz))

    ! Save Tn to 3D array for use with tavg_sgs
    Tn_all(:,:,jz) = Tn(:,:)
#ifdef ENABLE_CUDA
#ifndef PPDYN_TN
    end if
#endif
#endif
end do

! Share new data between overlapping nodes
#ifdef PPMPI
call mpi_sync_real_array( F_LM, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( F_MM, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( Tn_all, 0, MPI_SYNC_DOWNUP )
#ifdef PPDYN_TN
call mpi_sync_real_array( F_ee2, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( F_deedt2, 0, MPI_SYNC_DOWNUP )
call mpi_sync_real_array( ee_past, 0, MPI_SYNC_DOWNUP )
#endif
#endif

#ifdef PPLVLSET
call level_set_Cs_lag_dyn ()
#endif

! Reset variable for use during next set of cs_count timesteps
if( use_cfl_dt ) lagran_dt = 0._rprec

end subroutine lagrange_Ssim

#ifdef ENABLE_CUDA
!*******************************************************************************
logical function lagrange_ssim_cuda_enabled()
!*******************************************************************************
implicit none

lagrange_ssim_cuda_enabled = .true.

end function lagrange_ssim_cuda_enabled

!*******************************************************************************
subroutine lagrange_ssim_cuda_sync(where)
!*******************************************************************************
use cudafor
implicit none

character(len=*), intent(in) :: where
integer :: istat

istat = cudaDeviceSynchronize()
if (istat /= 0) then
    print *, 'lagrange_Ssim CUDA sync failure at ', trim(where), ': ', istat
    stop
end if
istat = cudaGetLastError()
if (istat /= 0) then
    print *, 'lagrange_Ssim CUDA kernel failure at ', trim(where), ': ', istat
    stop
end if

end subroutine lagrange_ssim_cuda_sync
#endif
