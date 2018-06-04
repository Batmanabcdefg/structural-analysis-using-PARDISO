
! (c)2018, Jakub Mikula
! mikula_jakub@centrum.sk

! PURPOSE:      index mesh analysis of a regular mesh on a rectangular domain
!               the topology (connectivity) matrix 
!
! SIZE of JA:
! 2*((na-2)*(nb-2)*9*2 + 2*(na-2)*6*2 + 2*(nb-2)*6*2 + 4*4*2)
! 
! ^       ^                  ^             ^           ^-- 4 corner nodes
! |       |                  |-------------|-- number of side nodes
! |       |                 
! |       |-- number of inner nodes
! |       
! |-- 2 dof in each node 
!
! the second number is the number of nearest neighbours + itself 
! the third number is 2 dofs per node 
! 
! ----------------------------------------------------------

subroutine mesh_analysis_CSR_2D(JA,IA,TFEM_CSR,K_diag,na,nb)

! DECLARE VARIABLES
! -----------------
INTEGER, intent(IN) :: na
INTEGER, intent(IN) :: nb

INTEGER, intent(OUT) :: JA(2*((na-2)*(nb-2)*9*2 + 2*(na-2)*6*2 + 2*(nb-2)*6*2 + 4*4*2))

INTEGER, intent(OUT) :: IA(na*nb*2+1)
INTEGER, intent(OUT) :: TFEM_CSR((na-1)*(nb-1),8,8)
INTEGER, intent(OUT) :: K_diag(2*na*nb)

INTEGER :: kkk, swch, pos, info, info_element
INTEGER :: perm_element(32)
INTEGER :: perm_dof1(32), perm_dof2(32)
INTEGER :: indx_i(32), indx_j(32)

INTEGER :: TFEM(4,8)
INTEGER :: x_e(4,8)
INTEGER :: y_e(4,8)

INTEGER :: k
INTEGER :: RR(14)
INTEGER :: RR2(18)
INTEGER :: RR2_excl(14)

INTEGER :: i1,i2
INTEGER :: dof, nxdof

INTEGER :: a_element, b_node
INTEGER :: element

! INITIALIZE VARIABLES
! --------------------
kkk     = 1
swch    = 1
pos     = 2
info    = 0

IA(1)   = 1

! degress of freedom numbering :
! o------o------o
! |13,14 |15,16 |17,18
! |      |      |
! o------o------o
! |7,8   |9,10  |11,12
! |      |      |
! o------o------o
! 1,2  3,4  5,6
!
! elements numbering:
! o------o------o
! |      |      | 
! | (3)  |  (4) |
! o------o------o
! |      |      |
! | (1)  |  (2) |
! o------o------o

! Topology matrix for the four elements in the COO format
! see figure above
TFEM(1,:) = (/1,2,3,4,9,10,7,8/)
TFEM(2,:) = (/3,4,5,6,11,12,9,10/)
TFEM(3,:) = (/7,8,9,10,15,16,13,14/)
TFEM(4,:) = (/9,10,11,12,17,18,15,16/)

! Unit element coordinates of the four elements
! see figure above
x_e(1,:) = (/-1,-1,0,0,0,0,-1,-1/)
x_e(2,:) = (/0,0,1,1,1,1,0,0/)
x_e(3,:) = (/-1,-1,0,0,0,0,-1,-1/)
x_e(4,:) = (/0,0,1,1,1,1,0,0/)

y_e(1,:) = (/-1,-1,-1,-1,0,0,0,0/)
y_e(2,:) = (/-1,-1,-1,-1,0,0,0,0/)
y_e(3,:) = (/0,0,0,0,1,1,1,1/)
y_e(4,:) = (/0,0,0,0,1,1,1,1/)

perm_element(:) = 0
perm_dof1(:) = 0
perm_dof2(:) = 0
indx_i(:) = 0
indx_j(:) = 0


! Fill in arrays 'perm_element', 'perm_dof2', 'indx_i', and 'indx_j'
! these arrays represent the permutation of a property as the conectivity of dof goes in the ascending order 
! for instance 11, 12, 13, 14, ... , 21, 22, 23, 24, ... etc. for the inner node of the four elements
! ------------------------------------------------------------------
k = 1
DO,i=1,18,1     !total number of degrees of freedom for the four elements
 DO,i1=1,4,1
 DO,i2=1,8,1
  IF (TFEM(i1,i2).EQ.i) THEN
    perm_element(k) = i1
    perm_dof2(k) = i2
    indx_i(k) = x_e(i1,i2)
    indx_j(k) = y_e(i1,i2)
    k = k + 1
  ENDIF
 ENDDO
 ENDDO
ENDDO !i

DO,i=1,32,1
 IF (perm_element(i) .EQ. 1) perm_dof1(i) = 5
 IF (perm_element(i) .EQ. 2) perm_dof1(i) = 7
 IF (perm_element(i) .EQ. 3) perm_dof1(i) = 3
 IF (perm_element(i) .EQ. 4) perm_dof1(i) = 1
ENDDO !i

! Find at what kkk to skip adding kkk, save into RR
k = 1
DO,i=1,32,1
 IF (perm_element(i) .NE. perm_element(i+1) .AND. perm_element(i) .LT. perm_element(i+1)) THEN
         RR(k) = i
         k = k + 1
 ENDIF
ENDDO !i

i=1
DO,dof=1,31,1
 IF (perm_element(dof) .NE. perm_element(dof+1)) THEN
  IF (perm_element(dof+1) >= perm_element(dof)) THEN
   IF (indx_i(dof) .EQ. indx_i(dof+1) .OR. indx_j(dof) .EQ. indx_j(dof+1)) THEN
        RR2_excl(i) = dof+1
        i = i + 1
   ENDIF
  ENDIF
 ENDIF
ENDDO !dof

i2 = 1
DO,i=1,32,1
 k = 0
 DO,i1=1,14
  IF (RR2_excl(i1) .EQ. i) THEN 
          k = 1
  ENDIF
 ENDDO !i1
 IF (k.EQ.0) THEN 
   RR2(i2) = i
   i2 = i2 + 1
 ENDIF
ENDDO !i

! ------------------------------------------------------------------------------
DO,j=1,nb,1
DO,i=1,na,1

 DO,nxdof=1,2,1
  DO,dof=1,32,1

   DO,k=1,18
        IF (RR2(k) == dof) THEN
         info = 0        
        ENDIF
   ENDDO !k

   info_element = 0
   SELECT CASE (perm_element(dof))
        CASE (1)
                IF (i>1 .AND. j>1) THEN
                 element = (j-2)*(na-1)+i-1
                 info = 1
                 info_element = 1
                ENDIF
        CASE (2)
                IF (i<na .AND. j>1) THEN
                 element = (j-2)*(na-1)+i
                 info = 1
                 info_element = 1
                ENDIF
        CASE (3)
                IF (i>1 .AND. j<nb) THEN
                 element = (j-1)*(na-1)+i-1
                 info = 1
                 info_element = 1
                ENDIF
        CASE (4)
                IF (i<na .AND. j<nb) THEN
                 element = (j-1)*(na-1)+i
                 info = 1
                 info_element = 1
                ENDIF
        CASE DEFAULT
   END SELECT

! ASSEMBLY THE csr TOPOLOGY MATRIX
! --------------------------------
   IF (info_element.EQ.1) TFEM_CSR(element,perm_dof1(dof)+nxdof-1,perm_dof2(dof)) = kkk

! FIND POSITION OF EACH NODE ON THE DIAGONAL OF THE STIFFNESS MATRIX
   IF (nxdof.EQ.1 .AND. dof.EQ.13) K_diag(2*(na*(j-1)+i)-1) = kkk
   IF (nxdof.EQ.2 .AND. dof.EQ.17) K_diag(2*(na*(j-1)+i)) = kkk

        k = 0
        DO,i1=1,14,1
         IF (RR(i1).EQ.dof) k = 1
        ENDDO

        IF (k.EQ.0 .AND. info.EQ.1) THEN
         JA(kkk) = 2*((j-1+indx_j(dof))*na+i+indx_i(dof))-1*swch
         IF (swch.EQ.0) THEN
                 swch = 1
         ELSE 
                 swch = 0
         ENDIF

        kkk = kkk + 1
        ENDIF

  ENDDO !dof

  IA(pos) = kkk
  pos = pos + 1
 ENDDO !nxdof

ENDDO !na
ENDDO !nb


end subroutine
