! (c) 2018, Jakub Mikula
! @ IHPC (Institute of High Performance Computing) A*STAR, Singapore
! ME (Faculty of Mechanical Engineering) NUS, Singapore
!
! CONTACT: mikula_jakub@centrum.sk
!
! LAST UPDATE: 20.08.2018

subroutine solve_PARDISO(pt,A,IA,JA,ndof,sizeA,F,x,phaseP,mtype)
USE omp_lib
INTEGER, intent(IN) :: ndof
INTEGER, intent(IN) :: sizeA
INTEGER :: i

INTEGER, intent(IN) :: phaseP(2)

DOUBLE PRECISION, intent(IN)   :: A(sizeA)
DOUBLE PRECISION, intent(IN)   :: F(ndof)
INTEGER, intent(IN)            :: JA(sizeA)
INTEGER, intent(IN)            :: IA(ndof+1)
DOUBLE PRECISION, intent(OUT)   :: x(ndof)

! PARDISO parameters
INTEGER*8                       :: pt(64)
INTEGER                         :: phase, error, msglvl
INTEGER                         :: iparm(64)
DOUBLE PRECISION                :: dparm(64)
INTEGER                         :: idum, solver
DOUBLE PRECISION                :: ddum
INTEGER                         :: nrhs, maxfct, mnum

INTEGER                         :: ogmt

INTEGER, intent(IN)             :: mtype

msglvl = 1
error = 0
solver = 1
maxfct = 1
mnum = 1
nrhs = 1

! follow 'Parallel Sparse Direct Solver PARDISO - User Guide Version 6.1.0'
iparm(1) = 0
iparm(2) = 3
iparm(3) = OMP_GET_MAX_THREADS()
!iparm(60) = 2
!iparm(4) = 31

! ----------------------------------------------------------------------------------
! Choose  phase of PARDISO solver

print*,size(F)
print*,size(x)
print*,size(IA)
print*,size(JA)

print*,ndof

IF (phaseP(1).EQ.-1) THEN
phase = -1 !termination and release of memory
print*,'Entered phase 11 (symbolic analysis)'
        CALL PARDISO(pt, maxfct, mnum, mtype, phase, ndof, A, IA, JA, idum,&
        nrhs, iparm, msglvl, F, x, error)
        IF (error.ne.0) THEN
                WRITE(*,*) 'The following PARDISO error detected (phase=11): ',error
                STOP 'algorithm terminated'
        ENDIF
print*,'Exited phase 11 (symbolic analysis)'
ENDIF

IF (phaseP(1).EQ.1) THEN
phase = 11 !symbolic analysis
print*,'Entered phase 11 (symbolic analysis)'
        CALL PARDISO(pt, maxfct, mnum, mtype, phase, ndof, A, IA, JA, idum,&
        nrhs, iparm, msglvl, F, x, error)
        IF (error.ne.0) THEN
                WRITE(*,*) 'The following PARDISO error detected (phase=11): ',error
                STOP 'algorithm terminated'
        ENDIF
print*,'Exited phase 11 (symbolic analysis)'
ENDIF

IF (phaseP(1)<= 2 .AND. phaseP(2)>=2) THEN
phase = 22 !Gauss elimination method
print*,'Entered phase 22 (Gauss elimination method)'
        CALL PARDISO(pt, maxfct, mnum, mtype, phase, ndof, A, IA, JA, idum,&
        nrhs, iparm, msglvl, F, x, error)
        IF (error.ne.0) THEN
                WRITE(*,*) 'The following PARDISO error detected (phase=22): ',error
                STOP 'algorithm terminated'
        ENDIF
print*,'Exited phase 22 (Gauss elimination method)'
ENDIF

IF (phaseP(1)<=3 .AND. phaseP(2)>=3) THEN
phase = 33 !backward substitution
print*,'Entered phase 33 (backward substitution)'
        CALL PARDISO(pt, maxfct, mnum, mtype, phase, ndof, A, IA, JA, idum,&
        nrhs, iparm, msglvl, F, x, error)
        IF (error.ne.0) THEN
                WRITE(*,*) 'The following PARDISO error detected (phase=33): ',error
                STOP 'algorithm terminated'
        ENDIF
print*,'Exited phase 33 (backward substitution)'
ENDIF

! Number of performed iterative refinement steps
!        OPEN(UNIT=1,FILE='iparm7.txt',POSITION='append')
!        WRITE(1,'(I5)') iparm(7)
!        CLOSE(1)


end subroutine
