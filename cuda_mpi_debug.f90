!*******************************************************************************
module cuda_mpi_debug
!*******************************************************************************
use mpi
use types, only : rprec
implicit none

save
private

public :: cuda_mpi_debug_init
public :: mpi_dbg_sendrecv_r, mpi_dbg_send_r, mpi_dbg_recv_r
public :: mpi_dbg_sendrecv_c

logical :: dbg_enabled = .false.
logical :: sync_enabled = .false.
logical :: dbg_initialized = .false.
integer :: dbg_count = 0
integer :: dbg_limit = 4000

contains

!*******************************************************************************
logical function cuda_mpi_env_set_not_zero(name)
!*******************************************************************************
character(*), intent(in) :: name
character(256) :: env
integer :: stat, lenv

call get_environment_variable(name, env, length=lenv, status=stat)
cuda_mpi_env_set_not_zero = (stat == 0 .and. lenv > 0 .and.                 &
    trim(env(1:lenv)) /= '0')

end function cuda_mpi_env_set_not_zero

!*******************************************************************************
subroutine cuda_mpi_debug_init(context, nproc_in, rank_in, coord_in, up_in, down_in, nz_in)
!*******************************************************************************
character(*), intent(in) :: context
integer, intent(in) :: nproc_in, rank_in, coord_in, up_in, down_in, nz_in
integer :: ierr_local
integer :: world_rank, world_size

if (dbg_initialized) return
dbg_initialized = .true.

dbg_enabled = cuda_mpi_env_set_not_zero('LESGO_MPI_CUDA_DEBUG')
sync_enabled = dbg_enabled

if (cuda_mpi_env_set_not_zero('LESGO_MPI_CUDA_SYNC')) sync_enabled = .true.

if (.not. dbg_enabled .and. .not. sync_enabled) return

call mpi_comm_rank(MPI_COMM_WORLD, world_rank, ierr_local)
call mpi_comm_size(MPI_COMM_WORLD, world_size, ierr_local)

if (.not. dbg_enabled) then
    write(*,*) 'CUDA_MPI_SYNC init context=', trim(context),                  &
        ' world_rank=', world_rank, ' world_size=', world_size,               &
        ' rank=', rank_in, ' coord=', coord_in, ' nproc=', nproc_in,          &
        ' up=', up_in, ' down=', down_in, ' nz=', nz_in
    call print_env('MPICH_GPU_SUPPORT_ENABLED')
    call print_env('CUDA_VISIBLE_DEVICES')
    return
endif

write(*,*) 'CUDA_MPI_DEBUG init context=', trim(context),                      &
    ' world_rank=', world_rank, ' world_size=', world_size,                    &
    ' rank=', rank_in, ' coord=', coord_in, ' nproc=', nproc_in,               &
    ' up=', up_in, ' down=', down_in, ' nz=', nz_in

call print_env('SLURM_NODEID')
call print_env('SLURM_LOCALID')
call print_env('SLURM_PROCID')
call print_env('SLURM_NTASKS')
call print_env('SLURM_NODELIST')
call print_env('CUDA_VISIBLE_DEVICES')
call print_env('MPICH_GPU_SUPPORT_ENABLED')
call print_env('MPICH_GPU_IPC_ENABLED')
call print_env('FI_PROVIDER')
call print_env('LD_LIBRARY_PATH')

write(*,*) 'CUDA_MPI_DEBUG legacy CUDA path removed'

end subroutine cuda_mpi_debug_init

!*******************************************************************************
subroutine print_env(name)
!*******************************************************************************
character(*), intent(in) :: name
character(1024) :: value
integer :: stat, lenv

call get_environment_variable(name, value, length=lenv, status=stat)
if (stat == 0) then
    if (lenv > 0) then
        write(*,*) 'CUDA_MPI_DEBUG env ', trim(name), '=', trim(value(1:lenv))
    else
        write(*,*) 'CUDA_MPI_DEBUG env ', trim(name), '='
    endif
else
    write(*,*) 'CUDA_MPI_DEBUG env ', trim(name), '=<unset>'
endif

end subroutine print_env

!*******************************************************************************
logical function should_print()
!*******************************************************************************
should_print = dbg_enabled .and. dbg_count < dbg_limit
if (should_print) dbg_count = dbg_count + 1
end function should_print

!*******************************************************************************
logical function should_sync_label(label)
!*******************************************************************************
character(*), intent(in) :: label

should_sync_label = sync_enabled

end function should_sync_label

!*******************************************************************************
subroutine cuda_pre(label, verbose)
!*******************************************************************************
character(*), intent(in) :: label
logical, intent(in) :: verbose
integer :: ierr_local, world_rank

if (.not. sync_enabled) return
if (verbose) call mpi_comm_rank(MPI_COMM_WORLD, world_rank, ierr_local)
if (verbose) write(*,*) 'CUDA_MPI_DEBUG PRE ', trim(label),                   &
    ' world_rank=', world_rank, ' no CUDA'

end subroutine cuda_pre

!*******************************************************************************
subroutine cuda_post(label, mpi_ierr, verbose)
!*******************************************************************************
character(*), intent(in) :: label
integer, intent(in) :: mpi_ierr
logical, intent(in) :: verbose
integer :: ierr_local, world_rank

if (.not. sync_enabled) return
if (verbose) call mpi_comm_rank(MPI_COMM_WORLD, world_rank, ierr_local)
if (verbose .or. mpi_ierr /= 0) write(*,*) 'CUDA_MPI_DEBUG POST ',            &
    trim(label), ' world_rank=', world_rank, ' mpi_ierr=', mpi_ierr, ' no CUDA'

end subroutine cuda_post

!*******************************************************************************
subroutine probe_real(label, buf)
!*******************************************************************************
character(*), intent(in) :: label
real(rprec), dimension(*) :: buf

if (.not. dbg_enabled) return

end subroutine probe_real

!*******************************************************************************
subroutine probe_complex(label, buf)
!*******************************************************************************
character(*), intent(in) :: label
complex(rprec), dimension(*) :: buf

if (.not. dbg_enabled) return

end subroutine probe_complex

!*******************************************************************************
subroutine mpi_dbg_sendrecv_r(sendbuf, sendcount, sendtype, dest, sendtag,      &
    recvbuf, recvcount, recvtype, source, recvtag, comm, status, ierr, label)
!*******************************************************************************
real(rprec), dimension(*) :: sendbuf
real(rprec), dimension(*) :: recvbuf
integer, intent(in) :: sendcount, sendtype, dest, sendtag
integer, intent(in) :: recvcount, recvtype, source, recvtag, comm
integer, intent(out) :: status(MPI_STATUS_SIZE)
integer, intent(out) :: ierr
character(*), intent(in) :: label
logical :: do_print, do_sync

do_print = should_print()
do_sync = should_sync_label(label)

if (do_print) then
    write(*,*) 'CUDA_MPI_DEBUG MPI_SENDRECV_R ', trim(label),                 &
        ' sendcount=', sendcount, ' dest=', dest, ' sendtag=', sendtag,       &
        ' recvcount=', recvcount, ' source=', source, ' recvtag=', recvtag
    call probe_real(trim(label)//' sendbuf', sendbuf)
    call probe_real(trim(label)//' recvbuf', recvbuf)
    call cuda_pre(trim(label)//' sendrecv', .true.)
else if (do_sync) then
    call cuda_pre(trim(label)//' sendrecv', .false.)
endif

call mpi_sendrecv(sendbuf, sendcount, sendtype, dest, sendtag,                &
    recvbuf, recvcount, recvtype, source, recvtag, comm, status, ierr)

if (do_print) then
    call cuda_post(trim(label)//' sendrecv', ierr, .true.)
else if (do_sync) then
    call cuda_post(trim(label)//' sendrecv', ierr, .false.)
endif

end subroutine mpi_dbg_sendrecv_r

!*******************************************************************************
subroutine mpi_dbg_sendrecv_c(sendbuf, sendcount, sendtype, dest, sendtag,      &
    recvbuf, recvcount, recvtype, source, recvtag, comm, status, ierr, label)
!*******************************************************************************
complex(rprec), dimension(*) :: sendbuf
complex(rprec), dimension(*) :: recvbuf
integer, intent(in) :: sendcount, sendtype, dest, sendtag
integer, intent(in) :: recvcount, recvtype, source, recvtag, comm
integer, intent(out) :: status(MPI_STATUS_SIZE)
integer, intent(out) :: ierr
character(*), intent(in) :: label
logical :: do_print, do_sync

do_print = should_print()
do_sync = should_sync_label(label)

if (do_print) then
    write(*,*) 'CUDA_MPI_DEBUG MPI_SENDRECV_C ', trim(label),                 &
        ' sendcount=', sendcount, ' dest=', dest, ' sendtag=', sendtag,       &
        ' recvcount=', recvcount, ' source=', source, ' recvtag=', recvtag
    call probe_complex(trim(label)//' sendbuf', sendbuf)
    call probe_complex(trim(label)//' recvbuf', recvbuf)
    call cuda_pre(trim(label)//' sendrecv', .true.)
else if (do_sync) then
    call cuda_pre(trim(label)//' sendrecv', .false.)
endif

call mpi_sendrecv(sendbuf, sendcount, sendtype, dest, sendtag,                &
    recvbuf, recvcount, recvtype, source, recvtag, comm, status, ierr)

if (do_print) then
    call cuda_post(trim(label)//' sendrecv', ierr, .true.)
else if (do_sync) then
    call cuda_post(trim(label)//' sendrecv', ierr, .false.)
endif

end subroutine mpi_dbg_sendrecv_c

!*******************************************************************************
subroutine mpi_dbg_send_r(buf, count, datatype, dest, tag, comm, ierr, label)
!*******************************************************************************
real(rprec), dimension(*) :: buf
integer, intent(in) :: count, datatype, dest, tag, comm
integer, intent(out) :: ierr
character(*), intent(in) :: label
logical :: do_print, do_sync

do_print = should_print()
do_sync = should_sync_label(label)

if (do_print) then
    write(*,*) 'CUDA_MPI_DEBUG MPI_SEND_R ', trim(label),                     &
        ' count=', count, ' dest=', dest, ' tag=', tag
    call probe_real(trim(label)//' sendbuf', buf)
    call cuda_pre(trim(label)//' send', .true.)
else if (do_sync) then
    call cuda_pre(trim(label)//' send', .false.)
endif

call mpi_send(buf, count, datatype, dest, tag, comm, ierr)

if (do_print) then
    call cuda_post(trim(label)//' send', ierr, .true.)
else if (do_sync) then
    call cuda_post(trim(label)//' send', ierr, .false.)
endif

end subroutine mpi_dbg_send_r

!*******************************************************************************
subroutine mpi_dbg_recv_r(buf, count, datatype, source, tag, comm, status, ierr, label)
!*******************************************************************************
real(rprec), dimension(*) :: buf
integer, intent(in) :: count, datatype, source, tag, comm
integer, intent(out) :: status(MPI_STATUS_SIZE)
integer, intent(out) :: ierr
character(*), intent(in) :: label
logical :: do_print, do_sync

do_print = should_print()
do_sync = should_sync_label(label)

if (do_print) then
    write(*,*) 'CUDA_MPI_DEBUG MPI_RECV_R ', trim(label),                     &
        ' count=', count, ' source=', source, ' tag=', tag
    call probe_real(trim(label)//' recvbuf', buf)
    call cuda_pre(trim(label)//' recv', .true.)
else if (do_sync) then
    call cuda_pre(trim(label)//' recv', .false.)
endif

call mpi_recv(buf, count, datatype, source, tag, comm, status, ierr)

if (do_print) then
    call cuda_post(trim(label)//' recv', ierr, .true.)
else if (do_sync) then
    call cuda_post(trim(label)//' recv', ierr, .false.)
endif

end subroutine mpi_dbg_recv_r

end module cuda_mpi_debug

