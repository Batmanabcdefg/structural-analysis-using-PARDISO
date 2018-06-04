program main3D

INTEGER :: na,nb,nc
INTEGER, allocatable :: JA(:), IA(:)
INTEGER, allocatable :: TFEM_CSR(:,:,:)
INTEGER, allocatable :: K_diag(:)

na = 2
nb = 2
nc = 2


allocate(JA(3*((na-2)*(nb-2)*(nc-2)*27*3 + 2*(na-2)*(nc-2)*18*3 + &
        2*(nb-2)*(nc-2)*18*3 + 2*(na-2)*(nb-2)*18*3 + 8*8*3)))


allocate(IA(na*nb*nc*3+1))
allocate(TFEM_CSR((na-1)*(nb-1)*(nc-1),24,24))
allocate(K_diag(3*na*nb*nc))

CALL mesh_analysis_CSR_3D(JA,IA,TFEM_CSR,K_diag,na,nb,nc)
!print*,IA

end program
