!(c) 2018, Jakub Mikula
!
! compile the code:
! shared memory architectures: make
! distributed memory architectures: make mpi
!
! run the code:
!   shared memory architectures: ./kuf
!   distributed memory architectures: mpirun -np nproc ./kuf_mpi
!   where nproc is the actual number of nodes
!
!   use: export OMP_NUM_THREADS=<number of threads per MPI process>
!        setenv OMP_NUM_THREADS <number of threads per MPI process>
!

program MAIN
implicit NONE

#ifdef Parallel
        include 'mpif.h'
#endif
!include 'mkl_cluster_sparse_solver.h'
!include 'mkl_cluster_sparse_solver.f90'

integer, parameter              :: dp=kind(0.d0)
integer                         :: nproc, ierr

double precision, allocatable   :: A(:), b(:), x(:)
integer, allocatable            :: JA(:), IA(:)

! PARDISO parameters
integer*8                       :: pt(64)
integer                         :: mtype, phase, error, msglvl
integer                         :: iparm(64)
double precision                :: dparm(64)
integer                         :: idum, solver
double precision                :: ddum
integer                         :: nrhs, maxfct, mnum, n, i
integer                         :: ii1
integer                         :: myid

#ifdef Parallel
        ! initialize MPI
        CALL MPI_INIT(ierr)
        ! use mpi_init_thread instead of mpi_init; recommened by:
        ! https://software.intel.com/en-us/mkl-developer-reference-c-parallel-direct-sparse-solver-for-clusters-interface

!        CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, ii1, ierr)
        CALL MPI_COMM_SIZE( mpi_comm_world, nproc, ierr)
        CALL MPI_COMM_RANK( mpi_comm_world, myid, ierr)
#endif

write(*,*), 'parallel process'

! allocate matrix
ALLOCATE(A(20)) 
ALLOCATE(IA(9))
ALLOCATE(JA(20))
ALLOCATE(x(8))
ALLOCATE(b(8))

A = (/7.0_dp,1.0_dp,2.0_dp,7.0_dp,-4.0_dp,8.0_dp,2.0_dp,1.0_dp,5.0_dp,7.0_dp,&
9.0_dp,-4.0_dp,7.0_dp,3.0_dp,5.0_dp,17.0_dp,11.0_dp,-3.0_dp,2.0_dp,5.0_dp/)
IA = (/1,5,8,10,12,13,16,18,21/)
JA = (/1,3,6,7,2,3,5,3,8,4,7,2,3,6,8,2,7,3,7,8/)

b(:) = 1.0_dp
x(:) = 0.0_dp

do,i=1,64
 pt(i) = 0
enddo

if (myid.eq.0) then
        msglvl = 1 !only process 0 prints PARDISO information
else
        msglvl = 0
endif

error = 0
mtype = 11
solver = 1
maxfct = 1
mnum = 1
nrhs = 1
n = 8

! for idetails on iparm visit: 
! https://software.intel.com/en-us/mkl-developer-reference-c-cluster-sparse-solver-iparm-parameter
do,i=1,64,1
 iparm(i) = 0
 dparm(i) = 0.0_dp
enddo


dparm(1)        = 300
dparm(2)        = 1e-6
dparm(3)        = 5000
dparm(4)        = 10
dparm(5)        = 1e-2
dparm(6)        = 5e-3
dparm(7)        = 10
dparm(8)        = 500
dparm(9)        = 25

#ifdef Parallel
write(*,*), 'calling cluster_sparse_solver' 
!iparm(1) = 0 --> set all iparm values to default
        phase = 11
        CALL cluster_sparse_solver(pt, maxfct, mnum, mtype, phase, n, A, IA, JA, idum,&
        nrhs, iparm, msglvl, b, x, mpi_comm_world, error)
        if (error.ne.0) then
                write(*,*) 'The following PARDISO error detected (phase=11): ',error
        !        stop 'algorithm terminated'
        endif

        phase = 22
        CALL cluster_sparse_solver(pt, maxfct, mnum, mtype, phase, n, A, IA, JA, idum,&
        nrhs, iparm, msglvl, b, x, mpi_comm_world, error)
        if (error.ne.0) then
                write(*,*) 'The following PARDISO error detected (phase=22): ',error
        !        stop 'algorithm terminated'
        endif

        phase = 33
        CALL cluster_sparse_solver(pt, maxfct, mnum, mtype, phase, n, A, IA, JA, idum,&
        nrhs, iparm, msglvl, b, x, mpi_comm_world, error)
        if (error.ne.0) then
                write(*,*) 'The following PARDISO error detected (phase=33): ',error
        !        stop 'algorithm terminated'
        endif

#endif

#ifndef Parallel
        iparm(1) = 1
        iparm(3) = 8

        phase = 11
        CALL PARDISO(pt, maxfct, mnum, mtype, phase, n, A, IA, JA, idum,&
        nrhs, iparm, msglvl, b, x, error)
        if (error.ne.0) then
                write(*,*) 'The following PARDISO error detected (phase=11): ',error
        !        stop 'algorithm terminated'
        endif

        phase = 22
        CALL PARDISO(pt, maxfct, mnum, mtype, phase, n, A, IA, JA, idum,&
        nrhs, iparm, msglvl, b, x, error)
        if (error.ne.0) then
                write(*,*) 'The following PARDISO error detected (phase=22): ',error
        !        stop 'algorithm terminated'
        endif

        phase = 33
        CALL PARDISO(pt, maxfct, mnum, mtype, phase, n, A, IA, JA, idum,&
        nrhs, iparm, msglvl, b, x, error)
        if (error.ne.0) then
                write(*,*) 'The following PARDISO error detected (phase=33): ',error
        !        stop 'algorithm terminated'
        endif

#endif


print*,'x: ',x

call sleep(5)

end program
