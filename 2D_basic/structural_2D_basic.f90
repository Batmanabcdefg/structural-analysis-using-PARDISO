! (c) 2018, Jakub Mikula
! @ IHPC (Institute of High Performance Computing) A*STAR, Singapore
! ME (Faculty of Mechanical Engineering) NUS, Singapore
!
! CONTACT: mikula_jakub@centrum.sk
!
! PURPOSE:
! Solve structural analysis (Kx = F + F0) on a 2D rectangular domain (na x nb)
! K - global stiffness matrix
! x - displacement field (u,v)
! F - external nodal force
! F0 - nodal force due to eigen (plastic) strain
!
! The algorithm is divided into 6 phases:
! (1) - generate mesh in the COO and CSR format (used to assembly the global stiffness matrix)
!       OUTPUT: x_e, y_e, TFEM+COO, TFEM_CSR, K_diag, JA, IA
! (2) - assembly the global stiffness matrix 
!       OUTPUT: A
! (4) - apply displacement boundary conditions using the penalty method
!       OUTPUT: A (updated)
! (5) - apply external nodal force boundary conditions
!       OUTPUT: F
! (6) - solve the structural analysis Kx = F + F0
!       OUTPUT: x
!
! PARDISO DOCUMENTATION: 
! Olaf Schenk and Klaus Gartner, PARDISO, User Guide Version 6.1.0 (Updated January 06, 2018)

subroutine structural_2D_basic(pt,x,na,nb,LA,LB,C_3x3,phase,phaseP,&
                A,IA,JA,TFEM_CSR,TFEM_COO,K_diag,x_e,y_e,theta_e,&
                penalty,ijk,NtN,bc,fc,F_external,F_plastic,nodal_force,mtype)
USE omp_lib
implicit NONE

INTEGER, PARAMETER :: dp=kind(0.d0)     !double precision

INTEGER*8, intent(INOUT) ::             pt(64) !slver memory pointer


INTEGER, intent(IN) ::                  na, nb
REAL(dp), intent(IN) ::                 LA, LB

DOUBLE PRECISION, intent(OUT) ::                 x(2*na*nb)

INTEGER, intent(INOUT) ::               JA( ((4*3*2)*2+(6*3*2)*(na-2))*(nb-2)+((8*2)*2+(12*2)*(na-2))*2 ) 
INTEGER, intent(INOUT) ::               IA( 2*na*nb + 1 )
REAL(dp), intent(INOUT) ::              A( ((4*3*2)*2+(6*3*2)*(na-2))*(nb-2)+((8*2)*2+(12*2)*(na-2))*2 )
INTEGER, intent(INOUT) ::               TFEM_CSR((na-1)*(nb-1),8,8)
INTEGER, intent(INOUT) ::               K_diag(2*na*nb)

REAL(dp), intent(INOUT) ::              x_e(4,(na-1)*(nb-1))
REAL(dp), intent(INOUT) ::              y_e(4,(na-1)*(nb-1))
INTEGER, intent(INOUT) ::               TFEM_COO((na-1)*(nb-1),4)


REAL(dp), intent(IN) ::                  C_3x3(3,3)
INTEGER, intent(IN) ::                  phase(2)
INTEGER, intent(IN) ::                  phaseP(2)

!REAL(dp), intent(IN) ::                 theta(na,nb)
REAL(dp), intent(IN) ::                 theta_e((na-1)*(nb-1))

REAL(dp), intent(IN) ::                 penalty

INTEGER, intent(IN) ::                  ijk
INTEGER, intent(IN) ::                  NtN

INTEGER, intent(IN) ::                  bc, fc

REAL(dp), intent(INOUT) ::              F_external(2*na*nb)
REAL(dp), intent(INOUT) ::              F_plastic(2*na*nb)

REAL(dp), intent(IN) ::                 nodal_force
INTEGER, intent(IN) ::                  mtype

! ---
REAL(dp) :: k_el(8,8)
REAL(dp) :: C_11, C_12, C_44

INTEGER :: i,j,k

REAL(dp) :: F_x(na*nb)
REAL(dp) :: F_y(na*nb)

REAL(dp) :: al, bl !(half-length of the rectangular element)

REAL(dp) :: C_3x3_rot(3,3)
REAL(dp) :: strain_0(3,1,4)

REAL(dp) :: F_elem(8)
REAL(dp) :: F(2*na*nb)

INTEGER :: nele

! -----------------------------------------------------------------------
! Generate mesh
IF (phase(1) == 1) THEN
        CALL mesh__COO_2D(x_e,y_e,TFEM_COO,na,nb,LA,LB)
print*,'na'
print*,'nb'
        CALL mesh__analysis_CSR_2Dx2(JA,IA,TFEM_CSR,K_diag,na,nb)
ENDIF

! -----------------------------------------------------------------------
! Assembly the stiffness matrix
IF (phase(1)<=2 .AND. phase(2)>=2) THEN

! ATTENTION: Due to efficiency I do not use REDUCTION of A   
! For small  matrices  this may  bring a  random  noise  in  nproc grid points of 
! the overlapping domains as two  or  more processes may  be  writing in the same
! element of the stiffness matrix. See the section on parallelization in the main 
! documentation.

nele = (na-1)*(nb-1)

A(:) = 0.0_dp
        !$OMP PARALLEL DO PRIVATE(k_el,i,j)
        DO,k=1,nele !number of elements
         DO,i=1,8 !node1
          DO,j=1,8 !node2
                C_11 = C_3x3(1,1) !in the crystal coordinate system
                C_12 = C_3x3(1,2) !in the crystal coordiante system
                C_44 = C_3x3(3,3) !in the crystal coordinate system
                ! The tensor is rotated within the subroutine

                CALL stiffness4(k_el, x_e(:,k),y_e(:,k),C_11,C_12,4.0_dp*C_44,theta_e(k)) 
                        ! 4.0_dp (engineering strain)
                        ! due to erlier subroutines 
                A(TFEM_CSR(k,i,j)) = A(TFEM_CSR(k,i,j)) + k_el(i,j)
          ENDDO
         ENDDO
        ENDDO
        !$OMP END PARALLEL DO
ENDIF

! -----------------------------------------------------------------------
! Update displacement boundary conditions
IF (phase(1)<=4 .AND. phase(2)>=4) THEN
        CALL bc__CSR_2Dx2(A,SIZE(A),K_diag,penalty,na,nb,bc,ijk)
ENDIF

! -----------------------------------------------------------------------
! Update external force boundary conditions
IF (phase(1)<=5 .AND. phase(2)>=5) THEN
        CALL force__external_2Dx2(F_external,na,nb,nodal_force,fc,ijk,NtN)
ENDIF

! -----------------------------------------------------------------------
! Solve
IF (phase(1)<=6 .AND. phase(2)>=6) THEN
        F =  F_external + F_plastic
        CALL solve_PARDISO(pt,A,IA,JA,2*na*nb,SIZE(A),F,x,phaseP,mtype)
ENDIF

end subroutine structural_2D_basic
