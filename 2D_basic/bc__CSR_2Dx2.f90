! (c) 2018, Jakub Mikula
! @ IHPC (Institute of High Performance Computing) A*STAR, Singapore
! ME (Faculty of Mechanical Engineering) NUS, Singapore


subroutine bc__CSR_2Dx2(A,sizeA,K_diag,penalty,na,nb,opt,ijk)
implicit NONE

INTEGER, PARAMETER :: dp=kind(0.d0)     !double precision


INTEGER, intent(IN) :: na,nb
INTEGER, intent(IN) :: ijk

INTEGER, intent(IN) :: sizeA

REAL(dp), intent(INOUT) :: A(sizeA)
INTEGER, intent(IN) :: K_diag(2*na*nb)
REAL(dp), intent(IN) :: penalty
INTEGER, intent(IN) :: opt

INTEGER :: i,j

SELECT CASE (opt)

CASE(1)

!           >-----<
!           |     |
!           |     |
!       (1) >-----< (4)
!           ^     ^
!          (2)   (3)  

        A(1) = A(1) + penalty             !(1) 
        A(10) = A(10) + penalty             !(2)
        A(24*na-20) = A(24*na-20) + penalty   !(3)
        !A(24*na-29) = penalty   !(4)


CASE(2)

!         > |-----| <
!         > |     | <
!         > |     | <
!         > |-----| <
!           ^     ^
!          (2)    (3)  

        A(10) = A(10) + penalty             !(2)
        A(24*na-20) = A(24*na-20) + penalty   !(3)

        ! Constraint left side in x
        !$OMP PARALLEL DO
        DO,i=1,nb,1
        A(K_diag(2*(na*(i-1)+1)-1)) = A(K_diag(2*(na*(i-1)+1)-1)) + penalty !left constraint in x
        ENDDO
        !$OMP END PARALLEL DO

        ! Constrain right side in x
        !$OMP PARALLEL DO
        DO,i=1,nb,1
        A(K_diag(2*(na*i)-1)) = A(K_diag(2*(na*i)-1)) + penalty !right constraint in x

        ENDDO
        !$OMP END PARALLEL DO

CASE DEFAULT
        STOP 'bc not chosen correctly *** process terminated ***'
END SELECT




end subroutine
