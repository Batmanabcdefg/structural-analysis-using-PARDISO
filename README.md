# structural-analysis-using-PARDISO

##Algorithm
 - read input file (system size, $\Delta x$)


### Mesh Generation
- symbolic analysis (generate CSR topology matrix indexing A,IA,JA arrays); processed on one processor; processed at the beginning of the transient analysis only
- stiffness matrix (A) assembly; parallel (OpenMP or MPI), assemblied at each time step of the transient analysis

<p align="center">
    <img src="https://github.com/MikulaJakub/structural-analysis-using-PARDISO/blob/master/Figures/numbering_3D.png" width="650"/>
</p>


