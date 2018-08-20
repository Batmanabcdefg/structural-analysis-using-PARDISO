! node = (j-1)*na + i

! one degree of freedom in each node
! (node)
! *element*
!
! (7)----(8)----(9)
!  | *3*  |  *4* |    
!  |      |      |
! (4)----(5)----(6)
!  | *1*  |  *2* |
!  |      |      |
! (1)----(2)----(3)
!
! sequence: element - node
! 1-1				1
! 1-2, 2-1			2
! 2-2				3
! 1-4, 3-1			4
! 1-3, 2-4, 3-2, 4-1		5
! 2-3, 4-2			6
! 3-4    			7
! 3-3, 4-4			8
! 4-3 				9

! two degrees of freedom in each node
! (node)
! *element*
!
! (13)(14)-(15)(16)--(17)(18)
!   |         |         |
!   | *3*     | *4*     |
!   |         |         |
! (7)(8)----(9)(10)--(11)(12)
!   |         |         |
!   | *1*     |  *2*    |
!   |         |         |
! (1)(2)----(3)(4)----(5)(6)
!
! sequence: element - node - degree of freedom		index_i 	index_j
! 1-1-1				1-1			-1		-1
! 1-1-2				1-2			-1		-1
! 1-2-1, 2-1-1			1-3, 2,1		0,0		-1,-1
! 1-2-2, 2-1-2			1-4, 2,2		0,0		-1,-1
! 2-2-1  			2-3			1		-1
! 2-2-2  			2-4			1		-1
! 1-4-1, 3-1-1			1-7,3-1			-1,-1		0,0
! 1-4-2, 3-1-2			1-8,3-2			-1,-1		0,0
! 1-3-1, 2-4-1, 3-2-1, 4-1-1	1-5, 2-7, 3-3, 4-1      0,0,0,0		0,0,0,0
! 1-3-2, 2-4-2, 3-2-2, 4-1-2	1-6, 2-8, 3-4, 4-2      0,0,0,0		0,0,0,0
! 2-3-1, 4-2-1		 	2-5, 4-3                1,1		0,0
! 2-3-2, 4-2-2			2-6, 4-4                1,1		0,0
! 3-4-1				3-7                     -1		1
! 3-4-2				3-8                     -1		1
! 3-3-1, 4-4-1			3-5, 4-7                0,0		1,1
! 3-3-2, 4-4-2			3-6, 4-8                0,0		1,1
! 4-3-1				4-5                     1		1
! 4-3-2				4-6                     1		1


subroutine mesh__analysis_CSR_2Dx2(JA,IA,TFEM_CSR,K_diag,na,nb)
implicit NONE

INTEGER, intent(IN) :: na,nb
INTEGER, intent(OUT) :: JA( ((4*3*2)*2+(6*3*2)*(na-2))*(nb-2)+((8*2)*2+(12*2)*(na-2))*2 )

INTEGER, intent(OUT) :: IA(2*na*nb+1)
INTEGER, intent(OUT) :: TFEM_CSR((na-1)*(nb-1),8,8)
INTEGER, intent(OUT) :: K_diag(2*na*nb)

integer, dimension(32) :: perm_element, perm_dof, perm_dof1
integer, dimension(32) :: index_i, index_j
integer :: new_row
integer :: switch = 1 

INTEGER :: info
INTEGER :: info_element

INTEGER :: k2
INTEGER :: k,i,j,l

INTEGER :: element

perm_dof1      = (/5,5,5,7,5,7,7,7,5,3,5,3,5,7,3,1,5,7,3,1,7,1,7,1,3,3,3,1,3,1,1,1/) !odd dofs
perm_element   = (/1,1,1,2,1,2,2,2,1,3,1,3,1,2,3,4,1,2,3,4,2,4,2,4,3,3,3,4,3,4,4,4/)
perm_dof       = (/1,2,3,1,4,2,3,4,7,1,8,2,5,7,3,1,6,8,4,2,5,3,6,4,7,8,5,7,6,8,5,6/)
index_i        = (/-1,-1,0,0,0,0,1,1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,1,1,1,1,-1,-1,0,0,0,0,1,1/)
index_j        = (/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1/)

l = 1
IA(1) = 1
new_row = 2
DO,j=1,nb,1
DO,i=1,na,1

DO,k2=1,2,1 !two degrees of freedom
 DO,k=1,32,1

  IF (ANY((/1,2,3,5,7,8,9,11,13,17,21,23,25,26,27,29,31,32,33/)==k)) THEN
   info=0
  ENDIF
  info_element=0 !if 0 element does not exist, if 1 element does exist

 select case (perm_element(k)) 
  case(1)
  ! element 1
  IF (i-1>=1 .AND. j-1>=1) THEN
   info=1
   info_element=1
   element = (j-2)*(na-1)+i-1 !(1)
  ENDIF 

  case(2)
  ! element 2
  IF (i+1<=na .AND. j-1>=1) THEN
   info=1
   info_element=1
   element = (j-2)*(na-1)+i !(2)
  ENDIF

  case(3)
 ! element 3
  IF (i-1>=1 .AND. j+1<=nb) THEN
   info=1
   info_element=1
   element = (j-1)*(na-1)+i-1 !(3)
  ENDIF

  case(4)
 ! element 4  
  IF (i+1<=na .AND. j+1<=nb) THEN
   info=1
   info_element=1
   element = (j-1)*(na-1)+i !(4)
  ENDIF
  case default
 info=0
 end select 

        IF (info_element==1) TFEM_CSR(element,perm_dof1(k)+k2-1,perm_dof(k)) = l

        ! Find position of each node on the diagoanl of the stiffness matrix
        IF (k2==1 .AND. k==13) THEN
        K_diag(2*(na*(j-1)+i)-1) = l
        ENDIF
        IF (k2==2 .AND. k==17) THEN
        K_diag(2*(na*(j-1)+i)) = l
        ENDIF

  IF (.NOT.ANY((/3,5,9,11,13,14,15,17,18,19,21,23,27,29/)==k)) THEN 
  !print*,perm_element(k),perm_dof(k),A(l),i,j
   IF (info .EQ. 1) THEN

     !!$OMP ORDERED
     JA(l) = 2*(((j+index_j(k))-1)*na + (i+index_i(k)))-1*switch 
     !!$OMP END ORDERED
!     JA(l) = 2*(((j+index_j(k))-1)*na + (i+index_i(k)))
     IF (switch == 0) THEN
      switch = 1
     ELSE
      switch = 0
     ENDIF

     l=l+1 !this must be here
!     print*,perm_element(k),perm_dof(k),A(l-1),l
   ENDIF
  ENDIF

 ENDDO

! print*,'new row',new_row
 IA(new_row) = l
 new_row = new_row + 1
ENDDO



ENDDO !i
ENDDO !j








end subroutine mesh__analysis_CSR_2Dx2
