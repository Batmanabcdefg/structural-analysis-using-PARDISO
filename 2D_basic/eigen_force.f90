!(c) 2017, Jakub Mikula
! PURPOSE:
! Deterine the odal force from the eigen strain epsilon_0
!
! Energy: 
! 1/2(epsilon-epsilon_0)^T * C * (epsilon-epsilon_0)
! epsilon - total strain
! epsilon_0 - eigen strain
! 

subroutine eigen_force(F_elem,C_TENSOR,strain_0,component,al,bl)

! F_elem - eigen force for k-th element (output)
! k - element
! c_TENSOR - symmetric tensor of elastic constants
! epsilon_0 - eigen strain
! component - number of eigen strain components
! al,bl - element size (half length)

integer :: ii,jj
integer :: i,j
integer :: component

double precision :: strain_B(3,8)
double precision :: c_TENSOR(3,3)
double precision :: eigen_VECTOR(3)
double precision, intent(OUT) :: F_elem(8)
double precision :: N_x(4,4,4)
double precision :: N_y(4,4,4)
double precision :: Ng(4,4,4)
double precision :: wg(4),eta(4),xi(4)
double precision :: al,bl
double precision :: strain_0(3,component,4)

wg(1) = 0.34785484513745385737
wg(2) = 0.65214515486254614263
wg(3) = wg(2) 
wg(4) = wg(1) 

xi(1) = -0.86113631159405257522
xi(2) = -0.33998104358485626480
xi(3) = -xi(2) 
xi(4) = -xi(1)

eta = xi

DO,ii=1,4,1
DO,jj=1,4,1
 Ng(1,ii,jj) = 1.0d0/4.0d0*(1.0d0-xi(ii))*(1.0d0-eta(jj))!shape function - node 1
 Ng(2,ii,jj) = 1.0d0/4.0d0*(1.0d0+xi(ii))*(1.0d0-eta(jj))!shape function - node 2
 Ng(3,ii,jj) = 1.0d0/4.0d0*(1.0d0+xi(ii))*(1.0d0+eta(jj))!shape function - node 3
 Ng(4,ii,jj) = 1.0d0/4.0d0*(1.0d0-xi(ii))*(1.0d0+eta(jj))!shape function - node 4

 N_x(1,ii,jj) = -1.0d0/4.0d0*(1.0d0-eta(jj))
 N_x(2,ii,jj) = 1.0d0/4.0d0*(1.0d0-eta(jj))
 N_x(3,ii,jj) = 1.0d0/4.0d0*(1.0d0+eta(jj))
 N_x(4,ii,jj) = -1.0d0/4.0d0*(1.0d0+eta(jj))

 N_y(1,ii,jj) = -1.0d0/4.0d0*(1.0d0-xi(ii))
 N_y(2,ii,jj) = -1.0d0/4.0d0*(1.0d0+xi(ii))
 N_y(3,ii,jj) = 1.0d0/4.0d0*(1.0d0+xi(ii))
 N_y(4,ii,jj) = 1.0d0/4.0d0*(1.0d0-xi(ii))
ENDDO
ENDDO

 DO,ii=1,4 !gauss points
 DO,jj=1,4 !gauss points

 strain_B(1,1) = N_x(1,ii,jj)/al
 strain_B(1,2) = 0.0d0
 strain_B(1,3) = N_x(2,ii,jj)/al
 strain_B(1,4) = 0.0d0
 strain_B(1,5) = N_x(3,ii,jj)/al
 strain_B(1,6) = 0.0d0
 strain_B(1,7) = N_x(4,ii,jj)/al
 strain_B(1,8) = 0.0d0

 strain_B(2,1) = 0.0d0
 strain_B(2,2) = N_y(1,ii,jj)/bl
 strain_B(2,3) = 0.0d0
 strain_B(2,4) = N_y(2,ii,jj)/bl
 strain_B(2,5) = 0.0d0
 strain_B(2,6) = N_y(3,ii,jj)/bl
 strain_B(2,7) = 0.0d0
 strain_B(2,8) = N_y(4,ii,jj)/bl

 strain_B(3,1) = N_y(1,ii,jj)/bl
 strain_B(3,2) = N_x(1,ii,jj)/al
 strain_B(3,3) = N_y(2,ii,jj)/bl
 strain_B(3,4) = N_x(2,ii,jj)/al
 strain_B(3,5) = N_y(3,ii,jj)/bl
 strain_B(3,6) = N_x(3,ii,jj)/al
 strain_B(3,7) = N_y(4,ii,jj)/bl
 strain_B(3,8) = N_x(4,ii,jj)/al

! eigen_VECTOR(1) = dot_product((/Ng(1,ii,jj),Ng(2,ii,jj),Ng(3,ii,jj),Ng(4,ii,jj)/),&
!(/phi_1(TFEM(k,1)),phi_1(TFEM(k,2)),phi_1(TFEM(k,3)),phi_1(TFEM(k,4))/))*varepsilon_xx1 + &
!dot_product((/Ng(1,ii,jj),Ng(2,ii,jj),Ng(3,ii,jj),Ng(4,ii,jj)/),&
!(/phi_2(TFEM(k,1)),phi_2(TFEM(k,2)),phi_2(TFEM(k,3)),phi_2(TFEM(k,4))/))*varepsilon_xx2 

! eigen_VECTOR(2) = dot_product((/Ng(1,ii,jj),Ng(2,ii,jj),Ng(3,ii,jj),Ng(4,ii,jj)/),&
!(/phi_1(TFEM(k,1)),phi_1(TFEM(k,2)),phi_1(TFEM(k,3)),phi_1(TFEM(k,4))/))*varepsilon_yy1 + &
! dot_product((/Ng(1,ii,jj),Ng(2,ii,jj),Ng(3,ii,jj),Ng(4,ii,jj)/),&
!(/phi_2(TFEM(k,1)),phi_2(TFEM(k,2)),phi_2(TFEM(k,3)),phi_2(TFEM(k,4))/))*varepsilon_yy2 

! eigen_VECTOR(3) = dot_product((/Ng(1,ii,jj),Ng(2,ii,jj),Ng(3,ii,jj),Ng(4,ii,jj)/),&
!(/phi_1(TFEM(k,1)),phi_1(TFEM(k,2)),phi_1(TFEM(k,3)),phi_1(TFEM(k,4))/))*varepsilon_xy1 + &
!dot_product((/Ng(1,ii,jj),Ng(2,ii,jj),Ng(3,ii,jj),Ng(4,ii,jj)/),&
!(/phi_2(TFEM(k,1)),phi_2(TFEM(k,2)),phi_2(TFEM(k,3)),phi_2(TFEM(k,4))/))*varepsilon_xy2

eigen_VECTOR(:) = 0.0d0
DO,j=1,3,1 !size of eigen_VECTOR
DO,i=1,component,1
 eigen_VECTOR(j) = eigen_VECTOR(j) + &
 dot_product((/Ng(1,ii,jj),Ng(2,ii,jj),Ng(3,ii,jj),Ng(4,ii,jj)/),strain_0(j,i,:))
ENDDO !component
ENDDO

  F_elem = F_elem + 1.0d0*MATMUL(MATMUL(TRANSPOSE(strain_B),C_TENSOR),eigen_VECTOR)*wg(ii)*wg(jj)*al*bl
 ENDDO
 ENDDO
!
end subroutine
