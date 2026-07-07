#ifndef PPCONVEC_GPU
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

!*******************************************************************************
subroutine convec
!*******************************************************************************
!
! Computes the rotation convective term in physical space
!       c = - (u X vort)
! Uses 3/2-rule for dealiasing for more info see Canuto 1991 Spectral Methods,
! chapter 7
!
use types, only : rprec
use param, only : BOGUS, coord, lbc_mom, lbz, ld_big, nproc, nx, nx2, ny,    &
    ny2, nz, sgs, ubc_mom
use sim_param, only : u, v, w, dudy, dudz, dvdx, dvdz, dwdx, dwdy
use sim_param, only : RHSx, RHSy, RHSz
use fft, only : back, back_big, forw, forw_big, padd, unpadd

implicit none

integer :: jz
integer :: jz_min
integer :: jzLo, jzHi, jz_max  ! added for full channel capabilities

real(rprec), save, allocatable, dimension(:,:,:) :: cc_big,                    &
    u_big, v_big, w_big, vort1_big, vort2_big, vort3_big
logical, save :: arrays_allocated = .false.

real(rprec) :: const


! Boundary vorticity index used in the wall cross-product terms below.
! LES/SGS uses the first interior uvp slab at the lower wall; DNS/no-SGS keeps
! the boundary slab.  The upper-wall term uses nz-1 in both cases.
if (sgs) then
    jzLo = 2
    jzHi = nz-1
else
    jzLo = 1
    jzHi = nz-1
endif

if( .not. arrays_allocated ) then
   allocate( cc_big( ld_big,ny2,nz ) )
   allocate( u_big(ld_big, ny2, lbz:nz) )
   allocate( v_big(ld_big, ny2, lbz:nz) )
   allocate( w_big(ld_big, ny2, lbz:nz) )
   allocate( vort1_big( ld_big,ny2,nz ) )
   allocate( vort2_big( ld_big,ny2,nz ) )
   allocate( vort3_big( ld_big,ny2,nz ) )
   arrays_allocated = .true.
endif

! Recall dudz, and dvdz are on UVP node for k=1 only
! So du2 does not vary from arg2a to arg2b in 1st plane (k=1)

! Loop through horizontal slices
! MPI: u_big, v_big needed at jz = 0, w_big not needed though
! MPI: could get u{1,2}_big
const = 1._rprec/(nx*ny)
do jz = lbz, nz
    ! use RHSx,RHSy,RHSz for temp storage
    RHSx(:,:,jz)=const*u(:,:,jz)
    RHSy(:,:,jz)=const*v(:,:,jz)
    RHSz(:,:,jz)=const*w(:,:,jz)

    ! do forward fft on normal-size arrays
    call dfftw_execute_dft_r2c(forw, RHSx(:,:,jz), RHSx(:,:,jz))
    call dfftw_execute_dft_r2c(forw, RHSy(:,:,jz), RHSy(:,:,jz))
    call dfftw_execute_dft_r2c(forw, RHSz(:,:,jz), RHSz(:,:,jz))

    ! zero pad: padd takes care of the oddballs
    call padd(u_big(:,:,jz), RHSx(:,:,jz))
    call padd(v_big(:,:,jz), RHSy(:,:,jz))
    call padd(w_big(:,:,jz), RHSz(:,:,jz))

    ! Back to physical space
    call dfftw_execute_dft_c2r(back_big, u_big(:,:,jz), u_big(:,:,jz))
    call dfftw_execute_dft_c2r(back_big, v_big(:,:,jz), v_big(:,:,jz))
    call dfftw_execute_dft_c2r(back_big, w_big(:,:,jz), w_big(:,:,jz))
end do


! Do the same thing with the vorticity
do jz = 1, nz
    ! if dudz, dvdz are on u-nodes for jz=1, then we need a special
    ! definition of the vorticity in that case which also interpolates
    ! dwdx, dwdy to the u-node at jz=1
    if ( (coord == 0) .and. (jz == 1) ) then

        select case (lbc_mom)
        ! Stress free
        case (0)
            RHSx(:, :, 1) = 0._rprec
            RHSy(:, :, 1) = 0._rprec

        ! Wall (all cases >= 1)
        case (1:)
            ! Wall-node vorticity: interpolate dwdy to the first uvp node.
            RHSx(:, :, 1) = const * ( 0.5_rprec * (dwdy(:, :, 1) +             &
                dwdy(:, :, 2))  - dvdz(:, :, 1) )
            ! Wall-node vorticity: interpolate dwdx to the first uvp node.
            RHSy(:, :, 1) = const * ( dudz(:, :, 1) -                          &
                0.5_rprec * (dwdx(:, :, 1) + dwdx(:, :, 2)) )

        end select
  endif

  if ( (coord == nproc-1) .and. (jz == nz) ) then

     select case (ubc_mom)

     ! Stress free
     case (0)

         RHSx(:, :, nz) = 0._rprec
         RHSy(:, :, nz) = 0._rprec

      ! No-slip and wall model
      case (1:)

         ! RHSx = vort1 at uvp nz-1, stored in the w-node nz work slot.
         RHSx(:, :, nz) = const * ( 0.5_rprec * (dwdy(:, :, nz-1) +            &
            dwdy(:, :, nz)) - dvdz(:, :, nz-1) )
         ! RHSy = vort2 at uvp nz-1, stored in the w-node nz work slot.
         RHSy(:, :, nz) = const * ( dudz(:, :, nz-1) -                         &
            0.5_rprec * (dwdx(:, :, nz-1) + dwdx(:, :, nz)) )

      end select
   endif

    ! Boundary slabs with explicit wall vorticity formulas were handled above;
    ! do not overwrite them with the general interior expression.
    if (.not.(coord==0 .and. jz==1) .and. .not. (ubc_mom>0 .and.               &
        coord==nproc-1 .and. jz==nz)  ) then
        RHSx(:,:,jz)=const*(dwdy(:,:,jz)-dvdz(:,:,jz))
        RHSy(:,:,jz)=const*(dudz(:,:,jz)-dwdx(:,:,jz))
    end if

    RHSz(:,:,jz)=const*(dvdx(:,:,jz)-dudy(:,:,jz))

    ! do forward fft on normal-size arrays
    call dfftw_execute_dft_r2c(forw, RHSx(:,:,jz), RHSx(:,:,jz))
    call dfftw_execute_dft_r2c(forw, RHSy(:,:,jz), RHSy(:,:,jz))
    call dfftw_execute_dft_r2c(forw, RHSz(:,:,jz), RHSz(:,:,jz))
    call padd(vort1_big(:,:,jz), RHSx(:,:,jz))
    call padd(vort2_big(:,:,jz), RHSy(:,:,jz))
    call padd(vort3_big(:,:,jz), RHSz(:,:,jz))

    ! Back to physical space
    ! FFTW backward transforms are unnormalized; const is applied when the
    ! padded-grid products are assembled below.
    call dfftw_execute_dft_c2r(back_big, vort1_big(:,:,jz), vort1_big(:,:,jz))
    call dfftw_execute_dft_c2r(back_big, vort2_big(:,:,jz), vort2_big(:,:,jz))
    call dfftw_execute_dft_c2r(back_big, vort3_big(:,:,jz), vort3_big(:,:,jz))
end do

! RHSx
! redefinition of const
const=1._rprec/(nx2*ny2)

if (coord == 0) then
    ! the cc's contain the normalization factor for the upcoming fft's
    cc_big(:,:,1)=const*(v_big(:,:,1)*(-vort3_big(:,:,1))&
       +0.5_rprec*w_big(:,:,2)*(vort2_big(:,:,jzLo)))
    ! jzLo selects the LES/DNS lower-wall vorticity slab described above.
    ! The 0.5 factor interpolates w(:,:,2) to the first uvp node above the
    ! wall; changing it is a wall-discretization change that needs validation.
    jz_min = 2
else
    jz_min = 1
end if

if (coord == nproc-1 ) then  ! channel
    ! the cc's contain the normalization factor for the upcoming fft's
    cc_big(:,:,nz-1)=const*(v_big(:,:,nz-1)*(-vort3_big(:,:,nz-1))&
        +0.5_rprec*w_big(:,:,nz-1)*(vort2_big(:,:,jzHi)))
    ! jzHi selects the upper-wall vorticity slab described above.
    ! The 0.5 factor interpolates w(:,:,nz-1) to the uvp node below the wall;
    ! changing it is a wall-discretization change that needs validation.

    jz_max = nz-2
else
    jz_max = nz-1
end if

do jz = jz_min, jz_max    !nz-1   ! channel
    cc_big(:,:,jz)=const*(v_big(:,:,jz)*(-vort3_big(:,:,jz))&
        +0.5_rprec*(w_big(:,:,jz+1)*(vort2_big(:,:,jz+1))&
        +w_big(:,:,jz)*(vort2_big(:,:,jz))))
end do

! Loop through horizontal slices
do jz=1,nz-1
    call dfftw_execute_dft_r2c(forw_big, cc_big(:,:,jz),cc_big(:,:,jz))
    ! un-zero pad
    ! note: cc_big is going into RHSx
    call unpadd(RHSx(:,:,jz),cc_big(:,:,jz))
    ! Back to physical space
    call dfftw_execute_dft_c2r(back, RHSx(:,:,jz), RHSx(:,:,jz))
end do

! RHSy
! const remains the padded-grid inverse transform normalization.
if (coord == 0) then
    ! the cc's contain the normalization factor for the upcoming fft's
    cc_big(:,:,1)=const*(u_big(:,:,1)*(vort3_big(:,:,1))&
        +0.5_rprec*w_big(:,:,2)*(-vort1_big(:,:,jzLo)))
    ! jzLo selects the LES/DNS lower-wall vorticity slab described above.
    ! The 0.5 factor interpolates w(:,:,2) to the first uvp node above the
    ! wall; changing it is a wall-discretization change that needs validation.
    jz_min = 2
else
    jz_min = 1
end if

if (coord == nproc-1) then   ! channel
    ! the cc's contain the normalization factor for the upcoming fft's
    cc_big(:,:,nz-1)=const*(u_big(:,:,nz-1)*(vort3_big(:,:,nz-1))&
        +0.5_rprec*w_big(:,:,nz-1)*(-vort1_big(:,:,jzHi)))
    ! jzHi selects the upper-wall vorticity slab described above.
    !--the 0.5 * w(:,:,nz-1) is the interpolation of w to the uvp node at nz-1
    !  below the wall

    jz_max = nz-2
else
    jz_max = nz-1
end if

do jz = jz_min, jz_max  !nz - 1   ! channel
   cc_big(:,:,jz)=const*(u_big(:,:,jz)*(vort3_big(:,:,jz))&
        +0.5_rprec*(w_big(:,:,jz+1)*(-vort1_big(:,:,jz+1))&
        +w_big(:,:,jz)*(-vort1_big(:,:,jz))))
end do

do jz=1,nz-1
    call dfftw_execute_dft_r2c(forw_big, cc_big(:,:,jz), cc_big(:,:,jz))
    ! un-zero pad
    ! note: cc_big is going into RHSy
    call unpadd(RHSy(:,:,jz), cc_big(:,:,jz))

    ! Back to physical space
    call dfftw_execute_dft_c2r(back, RHSy(:,:,jz), RHSy(:,:,jz))
end do

! RHSz

if (coord == 0) then
    ! There is no convective acceleration of w at wall or at top.
    !--not really true at wall, so this is an approximation?
    !  perhaps its OK since we dont solve z-eqn (w-eqn) at wall (its a BC)
    !--wrong, we do solve z-eqn (w-eqn) at bottom wall --pj
    !--earlier comment is also wrong, it is true that RHSz = 0 at both walls and
    ! slip BC
    cc_big(:,:,1)=0._rprec
    !! ^must change for Couette flow ... ?
    jz_min = 2
else
    jz_min = 1
end if

if (coord == nproc-1) then     ! channel
    ! There is no convective acceleration of w at wall or at top.
    !--not really true at wall, so this is an approximation?
    !  perhaps its OK since we dont solve z-eqn (w-eqn) at wall (its a BC)
    !--but now we do solve z-eqn (w-eqn) at top wall --pj
    !--earlier comment is also wrong, it is true that RHSz = 0 at both walls and
    ! slip BC
    cc_big(:,:,nz)=0._rprec
    !! ^must change for Couette flow ... ?
    jz_max = nz-1
else
    jz_max = nz-1   !! or nz ?       ! channel
end if

!#ifdef PPMPI
!  if (coord == nproc-1) then
!    cc_big(:,:,nz)=0._rprec ! according to JDA paper p.242
!    jz_max = nz - 1
!  else
!    jz_max = nz
!  endif
!#else
!  cc_big(:,:,nz)=0._rprec ! according to JDA paper p.242
!  jz_max = nz - 1
!#endif

! channel
do jz = jz_min, jz_max    !nz - 1
    cc_big(:,:,jz) = const*0.5_rprec*(                                         &
        (u_big(:,:,jz)+u_big(:,:,jz-1))*(-vort2_big(:,:,jz))                   &
        +(v_big(:,:,jz)+v_big(:,:,jz-1))*(vort1_big(:,:,jz)))
end do

! Loop through horizontal slices
do jz=1,nz !nz - 1
    call dfftw_execute_dft_r2c(forw_big,cc_big(:,:,jz),cc_big(:,:,jz))

    ! un-zero pad
    ! note: cc_big is going into RHSz!!!!
    call unpadd(RHSz(:,:,jz),cc_big(:,:,jz))

    ! Back to physical space
    call dfftw_execute_dft_c2r(back,RHSz(:,:,jz),   RHSz(:,:,jz))
end do

#ifdef PPMPI
#ifdef PPSAFETYMODE
RHSx(:, :, 0) = BOGUS
RHSy(:, :, 0) = BOGUS
RHSz(: ,:, 0) = BOGUS
#endif
#endif

!--top level is not valid
#ifdef PPSAFETYMODE
RHSx(:, :, nz) = BOGUS
RHSy(:, :, nz) = BOGUS
if(coord<nproc-1) RHSz(:, :, nz) = BOGUS
#endif

end subroutine convec
#endif
