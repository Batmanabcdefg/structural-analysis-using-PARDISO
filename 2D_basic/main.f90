! (c) 2018, Jakub Mikula
! @ IHPC (Institute of High Performance Computing) A*STAR, Singapore
! ME (Faculty of Mechanical Engineering) NUS, Singapore
!
! CONTACT: mikula_jakub@centrum.sk
!
! A simple example: structural analysis
! SOLVING:
!
!           ^ F
!           |
!  |--------|
!  |        |
!  |        |
! >|--------|
!  ^        ^
! u=0,v=0   v=0
!
!
! ASSUMPTIONS:
!  > cubic elasticity (isotropic elasticity)
!  > plane strain
!

program main
use variables
implicit NONE

! Initialize variables
! Mesh and geometry:
na = 10
nb = 10
LA = 10.0d0
LB = 10.0d0

! Material properties:
C_3x3(:,:) = 0.0d0
C_3x3(1,1) = 10.8d0
C_3x3(1,2) = 6.2d0
C_3x3(2,1) = 6.2d0
C_3x3(2,2) = 10.8d0
C_3x3(3,3) = (C_3x3(1,1)-C_3x3(1,2))/2.0d0

! PARDISO
pt(:) = 0.0d0

! Allocate arrays
allocate(x(2*na*nb))
allocate(A(((4*3*2)*2+(6*3*2)*(na-2))*(nb-2)+((8*2)*2+(12*2)*(na-2))*2))
allocate(JA(SIZE(A)))
allocate(IA(2*na*nb+1))
allocate(TFEM_CSR((na-1)*(nb-1),8,8))
allocate(TFEM_COO((na-1)*(nb-1),4))
allocate(x_e(4,(na-1)*(nb-1)))
allocate(y_e(4,(na-1)*(nb-1)))
allocate(theta_e((na-1)*(nb-1)))

allocate(K_diag(2*na*nb))

allocate(F_external(2*na*nb))
allocate(F_plastic(2*na*nb))

! Initialize arrays
theta_e(:) = 0.0d0

CALL structural_2D_basic(pt,x,na,nb,LA,LB,C_3x3,(/1,6/),(/1,3/),&
A,IA,JA,TFEM_CSR,TFEM_COO,K_diag,x_e,y_e,theta_e,&
10000.0d0,1,1,1,1,F_external,F_plastic,10.0d0,11)


! Postprocessing
! Output the nodal displacements
open(unit=1,file='x.txt')
write(1,*) x
close(1)

end program 
