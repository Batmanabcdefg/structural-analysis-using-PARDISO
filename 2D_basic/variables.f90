! (c) 2018, Jakub Mikula
! @ IHPC (Institute of High Performance Computing) A*STAR, Singapore
! ME (Faculty of Mechanical Engineering) NUS, Singapore
!
! CONTACT: mikula_jakub@centrum.sk
!
module variables
INTEGER, PARAMETER :: dp=kind(0.d0) !double precision

INTEGER*8 :: pt(64)

INTEGER :: na, nb
REAL(dp) :: LA,LB
REAL(dp) :: C_3x3(3,3)

REAL(dp), allocatable :: x(:)
REAL(dp), allocatable :: A(:)
INTEGER, allocatable :: IA(:)
INTEGER, allocatable :: JA(:)
INTEGER, allocatable :: TFEM_CSR(:,:,:)
INTEGER, allocatable :: TFEM_COO(:,:)

REAL(dp), allocatable :: x_e(:,:)
REAL(dp), allocatable :: y_e(:,:)
!REAL(dp), allocatable :: theta(:,:)
REAL(dp), allocatable :: theta_e(:)

INTEGER, allocatable :: K_diag(:)

REAL(dp), allocatable :: F_external(:)
REAL(dp), allocatable :: F_plastic(:)

end module variables
