# structural-analysis-using-PARDISO

## Mesh Generation
- symbolic analysis (generate CSR topology matrix indexing A,IA,JA arrays); processed on one processor; processed at the beginning of the transient analysis only
- stiffness matrix (A) assembly; parallel (OpenMP or MPI), assemblied at each time step of the transient analysis
