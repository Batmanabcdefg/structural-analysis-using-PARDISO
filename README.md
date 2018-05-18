
# structural-analysis-using-PARDISO

## Algorithm
 - read input file (LA,LB,na,nb)
 - run symbolic mesh analysis to obtain: **IA**, **JA**, **TFEM_CSR**, **K_diag**
 - assembly the stiffness matrix **K**
 - provide the nodal force **F** and boundary conditions using the penalty method
 - solve **K** **u** = **F** using PARDISO
 - postprocess results (displacement field) **u**


### Mesh Generation
- symbolic analysis (generate CSR topology matrix indexing A,IA,JA arrays); processed on one processor; processed at the beginning of the transient analysis only
- stiffness matrix (A) assembly; parallel (OpenMP or MPI), assemblied at each time step of the transient analysis

<p align="center">
    <img src="https://github.com/MikulaJakub/structural-analysis-using-PARDISO/blob/master/Figures/numbering_3D.png" width="650"/>
</p>


