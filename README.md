
# structural-analysis-using-PARDISO

## Periodic Boundary conditions
Periodic boundary conditions are activated by initializing `periodic_constraint .TRUE.`. Then, specify the constraint. This is done by the vector `pphase(/s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9)/)`. The parameters `s(1)` to `s(9)` can be initialized with `.TRUE.` or `.FALSE.` depending on whether the constrains `s(i)` is active. `s(1)` - constrain x=0 in x; `s(2)` - constrain x=0 in y; `s(3)` -constrain x=0 in z; `s(4)` - constrain y=0 in x; etc. 

## Algorithm
 1) read input file (LA,LB,na,nb)
 2) run symbolic mesh analysis to obtain: **IA**, **JA**, **TFEM_CSR**, **K_diag**
 3) assembly the stiffness matrix **K** (*consider parallelization*)
 4) provide the nodal force **F** and boundary conditions using the penalty method
 5) solve **K** **u** = **F** using PARDISO (*consider parallelization*)
 6) postprocess results (displacement field) **u**


### 2) symbolic mesh analysis
It is assumed that the numbering of nodes, elements and degrees of freedom is such as shown in the figure below:

- 2D:

- 3D:
<p align="center">
    <img src="https://github.com/MikulaJakub/structural-analysis-using-PARDISO/blob/master/Figures/numbering_3D.png" width="650"/>
</p>

The red node is the node which connets all eight elements togehter. Each node contains three degrees of freedom u_x, u_y, and u_z. 

### 3) stiffness matrix assembly

