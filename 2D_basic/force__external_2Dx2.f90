! (c) 2018, Jakub Mikula
! @ IHPC (Institute of High Performance Computing) A*STAR, Singapore
! ME (Faculty of Mechanical Engineering) NUS, Singapore


subroutine force__external_2Dx2(F_external,na,nb,f,opt,ijk,NtN)
implicit NONE

INTEGER, PARAMETER :: dp=kind(0.d0)     !double precision

INTEGER, intent(IN) :: NtN

INTEGER, intent(IN) :: na,nb
INTEGER, intent(IN) :: ijk

REAL(dp), intent(INOUT) :: F_external(2*na*nb)

REAL(dp), intent(IN) :: f
INTEGER, intent(IN) :: opt

INTEGER :: k

INTEGER :: i

SELECT CASE (opt)

CASE(0)
      F_external = 0.0_dp

CASE(1)
      F_external(2*na*nb) = f

CASE(3)

      DO,i=1,nb,1
      F_external(2*(na*i)) = -0.05_dp*ijk/NtN !up
      F_external(2*(na*i)) = -2.0_dp*(-0.05)*ABS(ijk-NtN/2.0_dp)/NtN+(-0.05) !up and down
      ENDDO







CASE DEFAULT
      STOP 'fc not chosen correctly *** process terminated ***'
END SELECT






end subroutine force__external_2Dx2
