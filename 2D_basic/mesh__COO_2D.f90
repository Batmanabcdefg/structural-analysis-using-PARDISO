subroutine mesh__COO_2D(x_e,y_e,TFEM,ndof_a,ndof_b,La,Lb)
implicit NONE

INTEGER, PARAMETER :: dp=kind(0.d0)     !double precision


INTEGER, intent(IN) :: ndof_a,ndof_b
REAL(dp), intent(IN) :: La, Lb

REAL(dp), intent(OUT) :: x_e(4,(ndof_a-1)*(ndof_b-1))
REAL(dp), intent(OUT) :: y_e(4,(ndof_a-1)*(ndof_b-1))
INTEGER, intent(OUT) :: TFEM((ndof_a-1)*(ndof_b-1),4)

INTEGER :: nele, ndof
INTEGER :: i,j


nele = (ndof_a-1)*(ndof_b-1)
ndof = ndof_a*ndof_b
      
        ! Numbering of nodes:
        ! 0-----0-----0-----0-----0
        ! |(11) |(12) |(13) |(14) |(15)
        ! 0-----0-----0-----0-----0
        ! |(6)  |(7)  |(8)  |(9)  |(10)
        ! 0-----0-----0-----0-----0
        !  (1)   (2)   (3)   (4)   (!5)
        
         DO,j=1,ndof_b-1,1
         DO,i=1,ndof_a-1,1
          x_e(1,(ndof_a-1)*(j-1)+i) = La/(ndof_a-1)*(i-1)
          x_e(2,(ndof_a-1)*(j-1)+i) = La/(ndof_a-1)*i
          x_e(3,(ndof_a-1)*(j-1)+i) = La/(ndof_a-1)*i
          x_e(4,(ndof_a-1)*(j-1)+i) = La/(ndof_a-1)*(i-1)

          y_e(1,(ndof_a-1)*(j-1)+i) = Lb/(ndof_b-1)*(j-1)
          y_e(2,(ndof_a-1)*(j-1)+i) = Lb/(ndof_b-1)*(j-1)
          y_e(3,(ndof_a-1)*(j-1)+i) = Lb/(ndof_b-1)*j
          y_e(4,(ndof_a-1)*(j-1)+i) = Lb/(ndof_b-1)*j

          ! TFEM - one degree of freedom in each node
          TFEM((ndof_a-1)*(j-1)+i,1) = (j-1)*ndof_a + i 
          TFEM((ndof_a-1)*(j-1)+i,2) = (j-1)*ndof_a + i+1
          TFEM((ndof_a-1)*(j-1)+i,3) = j*ndof_a + i+1
          TFEM((ndof_a-1)*(j-1)+i,4) = j*ndof_a + i
         ENDDO
         ENDDO

end subroutine mesh__COO_2D
