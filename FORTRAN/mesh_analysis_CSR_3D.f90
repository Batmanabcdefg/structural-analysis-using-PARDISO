
! (c)2018, Jakub Mikula
! mikula_jakub@centrum.sk

! This project was supported through the SINGA schoolarship award
! 
! PURPOSE:      index mesh analysis of a regular mesh on a tetragonal domain
!               calculate the topology (connectivity) matrix
!
! SIZE fo JA:
! 3*((na-2)*(nb-2)*(nc-2)*27*3 + 3*(na-2)*18*3 + 3*(nb-2)*18*3 + 3*(nc-2)*18*3 + 12*8*3)
! 
! ^     ^                             ^               ^               ^           ^
! |     |                             |               |               |           |-- 12 corner nodes
! |     |                             |               |               |           
! |     |                             |---------------|---------------|-- number of side nodes
! |     |
! |     |-- number of inner nodes
! |
! |-- 3 dof in each node
!
! the second number os the number of nearest neighbours + itself
! the third number os the 3 dofs per node
!
! ----------------------------------------------------------
subroutine mesh_analysis_CSR_3D(JA,IA,TFEM_CSR,K_diag,na,nb,nc)

! DECLARE VARIABLES
! -----------------

INTEGER, intent(IN) :: na
INTEGER, intent(IN) :: nb
INTEGER, intent(IN) :: nc

INTEGER, intent(OUT) :: JA(3*((na-2)*(nb-2)*(nc-2)*27*3 + 2*(na-2)*(nc-2)*18*3 + &
        2*(nb-2)*(nc-2)*18*3 + 2*(na-2)*(nb-2)*18*3 + 8*8*3))
INTEGER, intent(OUT) :: IA(na*nb*nc*3+1)
INTEGER, intent(OUT) :: TFEM_CSR((na-1)*(nb-1)*(nc-1),24,24)
INTEGER, intent(OUT) :: K_diag(3*na*nb*nc)

INTEGER :: i,j,k,l !loop integers

INTEGER :: kkk, swch, pos, info, info_element, indx
INTEGER :: perm_element(192)
INTEGER :: perm_dof1(192), perm_dof2(192)
INTEGER :: indx_i(192), indx_j(192), indx_k(192)

INTEGER :: TFEM(8,24)
INTEGER :: x_e(8,24), y_e(8,24), z_e(8,24)

INTEGER :: RR(111)
INTEGER :: RR2(81)
INTEGER :: RR2_excl(111)

INTEGER :: dof, nxdof
INTEGER :: element
INTEGER :: i1,i2

! INITIALIZE VARIABLES
! --------------------
kkk     = 1
swch    = 2
pos     = 2
info    = 0
indx    = 2

IA(1)   = 1

perm_element(:) = 0
perm_dof1(:) = 0
perm_dof2(:) = 0
indx_i(:) = 0
indx_j(:) = 0
indx_k(:) = 0

TFEM(1,:) = (/1,2,3,4,5,6,13,14,15,10,11,12,28,29,30,31,32,33,40,41,42,37,38,39/)
TFEM(2,:) = (/4,5,6,7,8,9,16,17,18,13,14,15,31,32,33,34,35,36,43,44,45,40,41,42/)
TFEM(3,:) = (/10,11,12,13,14,15,22,23,24,19,20,21,37,38,39,40,41,42,49,50,51,46,47,48/)
TFEM(4,:) = (/13,14,15,16,17,18,25,26,27,22,23,24,40,41,42,43,44,45,52,53,54,49,50,51/)
TFEM(5,:) = (/28,29,30,31,32,33,40,41,42,37,38,39,55,56,57,58,59,60,67,68,69,64,65,66/)
TFEM(6,:) = (/31,32,33,34,35,36,43,44,45,40,41,42,58,59,60,61,62,63,70,71,72,67,68,69/)
TFEM(7,:) = (/37,38,39,40,41,42,49,50,51,46,47,48,64,65,66,67,68,69,76,77,78,73,74,75/)
TFEM(8,:) = (/40,41,42,43,44,45,52,53,54,49,50,51,67,68,69,70,71,72,79,80,81,76,77,78/)

x_e(1,:) = (/-1,-1,-1,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,-1,-1,-1/)
x_e(2,:) = (/0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0/)
x_e(3,:) = (/-1,-1,-1,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,-1,-1,-1/)
x_e(4,:) = (/0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0/)
x_e(5,:) = (/-1,-1,-1,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,-1,-1,-1/)
x_e(6,:) = (/0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0/)
x_e(7,:) = (/-1,-1,-1,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,-1,-1,-1/)
x_e(8,:) = (/0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0/)

y_e(1,:) = (/-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0/)
y_e(2,:) = (/-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0/)
y_e(3,:) = (/0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1/)
y_e(4,:) = (/0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1/)
y_e(5,:) = (/-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0/)
y_e(6,:) = (/-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0/)
y_e(7,:) = (/0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1/)
y_e(8,:) = (/0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1/)

z_e(1,:) = (/-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0/)
z_e(2,:) = (/-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0/)
z_e(3,:) = (/-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0/)
z_e(4,:) = (/-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0/)
z_e(5,:) = (/0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1/)
z_e(6,:) = (/0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1/)
z_e(7,:) = (/0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1/)
z_e(8,:) = (/0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1/)

l = 1
DO,i=1,81,1 !total number of degrees of freedom for the eight elements
 DO,i1=1,8,1
 DO,i2=1,24,1
  IF (TFEM(i1,i2).EQ.i) THEN
          perm_element(l) = i1
          perm_dof2(l) = i2
          indx_i(l) = x_e(i1,i2)
          indx_j(l) = y_e(i1,i2)
          indx_k(l) = z_e(i1,i2)
          l = l + 1
  ENDIF
 ENDDO !i2
 ENDDO !i1
ENDDO !i

DO,i=1,192,1
 IF (perm_element(i) .EQ. 1) perm_dof1(i) = 19
 IF (perm_element(i) .EQ. 2) perm_dof1(i) = 22
 IF (perm_element(i) .EQ. 3) perm_dof1(i) = 16
 IF (perm_element(i) .EQ. 4) perm_dof1(i) = 13
 IF (perm_element(i) .EQ. 5) perm_dof1(i) = 7
 IF (perm_element(i) .EQ. 6) perm_dof1(i) = 10
 IF (perm_element(i) .EQ. 7) perm_dof1(i) = 4
 IF (perm_element(i) .EQ. 8) perm_dof1(i) = 1
ENDDO !i

l = 1
DO,i=1,191,1
 IF (perm_element(i).NE.perm_element(i+1) .AND. perm_element(i).LT.perm_element(i+1)) THEN
         RR(l) = i
         l = l + 1
 ENDIF
ENDDO !i

l = 1
DO,dof=1,191,1
 IF (perm_element(dof) .NE. perm_element(dof+1)) THEN
  IF (perm_element(dof+1) >= perm_element(dof)) THEN
   IF (indx_i(dof).EQ.indx_i(dof+1) .OR. indx_j(dof).EQ.indx_j(dof+1) .OR. indx_k(dof).EQ.indx_k(dof+1)) THEN
           RR2_excl(l) = dof+1
           l = l + 1
   ENDIF
  ENDIF
 ENDIF
ENDDO !dof

i2 = 1
DO,i=1,192,1
 k = 0
 DO,i1=1,111
  IF (RR2_excl(i1).EQ.i) THEN
          k = 1
  ENDIF
 ENDDO !i1
 IF (k.EQ.0) THEN
         RR2(i2) = i
         i2 = i2 + 1
 ENDIF
ENDDO !i

! ------------------------------------------------------------------------------
DO,k=1,nc,1
DO,j=1,nb,1
DO,i=1,na,1

 DO,nxdof=1,3,1
 DO,dof=1,192,1

  DO,l=1,81,1
        IF (RR2(l) .EQ. dof) info = 0
  ENDDO !l

  info_element = 0
  SELECT CASE (perm_element(dof))
        CASE (1)
         IF (i>1 .AND. j>1 .AND. k>1) THEN
                element = (k-2)*(na-1)*(nb-1) + (j-2)*(na-1) + i-1
                info = 1
                info_element = 1
         ENDIF
        CASE (2)
         IF (i<na .AND. j>1 .AND. k>1) THEN
                element = (k-2)*(na-1)*(nb-1) + (j-2)*(na-1) + i
                info = 1
                info_element = 1
         ENDIF
        CASE (3)
         IF (i>1 .AND. j<nb .AND. k>1) THEN
                element = (k-2)*(na-1)*(nb-1) + (j-1)*(na-1) + i-1
                info = 1
                info_element = 1
         ENDIF
        CASE (4)
         IF (i<na .AND. j<nb .AND. k>1) THEN
                element = (k-2)*(na-1)*(nb-1) + (j-1)*(na-1) + i
                info = 1
                info_element = 1
         ENDIF
        CASE (5)
         IF (i>1 .AND. j>1 .AND. k<nc) THEN
                element = (k-1)*(na-1)*(nb-1) + (j-2)*(na-1) + i-1
                info = 1
                info_element = 1
         ENDIF
        CASE (6)
         IF (i<na .AND. j>1 .AND. k<nc) THEN
                element = (k-1)*(na-1)*(nb-1) + (j-2)*(na-1) + i
                info = 1
                info_element = 1
         ENDIF
        CASE (7)
         IF (i>1 .AND. j<nb .AND. k<nc) THEN
                element = (k-1)*(na-1)*(nb-1) + (j-1)*(na-1) + i-1
                info = 1
                info_element = 1
         ENDIF
        CASE (8)
         IF (i<na .AND. j<nb .AND. k<nc) THEN
                element = (k-1)*(na-1)*(nb-1) + (j-1)*(na-1) + i
                info = 1
                info_element = 1
         ENDIF
        CASE DEFAULT
 END SELECT

! ASSEMBLY THE csr TOPOLOGY MATRIX
 IF (info_element.EQ.1) TFEM_CSR(element,perm_dof1(dof)+nxdof-1,perm_dof2(dof)) = kkk

! FIND POSITION OF EACH NODE ON THE DIAGONAL OF THE STIFFNESS MATRIX
 IF (nxdof.EQ.1 .AND. dof.EQ.85) K_diag(3*(na*nb*(k-1)+na*(j-1)+i)-2) = kkk
 IF (nxdof.EQ.2 .AND. dof.EQ.93) K_diag(3*(na*nb*(k-1)+na*(j-1)+i)-1) = kkk
 IF (nxdof.EQ.3 .AND. dof.EQ.101) K_diag(3*(na*nb*(k-1)+na*(j-1)+i)) = kkk

  l = 0
  DO,i1=1,111,1
   IF (RR(i1).EQ.dof) l = 1
  ENDDO
  IF (l.EQ.0 .AND. info.EQ.1) THEN
          JA(kkk) = 3*((k-1+indx_k(dof))*na*nb + (j-1+indx_j(dof))*na + i+indx_i(dof)) -1*swch
          SELECT CASE (swch)
                CASE (0)
                        swch = 2
                CASE (1)
                        swch = 0
                CASE (2)
                        swch = 1
                CASE DEFAULT
          END SELECT
          kkk = kkk + 1
   ENDIF

 ENDDO !dof
  IA(pos) = kkk
  pos = pos + 1
 ENDDO !nxdof

ENDDO !i
ENDDO !j
ENDDO !k

print*,JA
end subroutine
