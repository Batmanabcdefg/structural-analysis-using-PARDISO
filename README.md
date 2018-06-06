
# structural-analysis-using-PARDISO

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

