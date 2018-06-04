program main

INTEGER :: na,nb
INTEGER, allocatable :: JA(:), IA(:)
INTEGER, allocatable :: TFEM_CSR(:,:,:)
INTEGER, allocatable :: K_diag(:)

na = 2
nb = 2

allocate(JA(2*((na-2)*(nb-2)*9*2 + 2*(na-2)*6*2 + 2*(nb-2)*6*2 + 4*4*2)))
allocate(IA(na*nb*2+1))
allocate(TFEM_CSR((na-1)*(nb-1),8,8))
allocate(K_diag(2*na*nb))

      CALL mesh_analysis_CSR_2D(JA,IA,TFEM_CSR,K_diag,na,nb) 

      print*,TFEM_CSR
end program
