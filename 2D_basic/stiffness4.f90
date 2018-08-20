
SUBROUTINE stiffness4(k_el, x_e,y_e,c_11_i,c_12_i,c_44_i,theta_e)
IMPLICIT NONE

!declare variables

    double precision :: theta_e
        
    integer :: i
    integer :: j
        
        double precision, intent(IN) :: C_11_i            !cubic elasticity 
        double precision, intent(IN) :: C_12_i            !cubic elasticity 
        double precision, intent(IN) :: C_44_i            !cubic elasticity 
        
        double precision :: C_11
        double precision :: C_12
        double precision :: C_13
        double precision :: C_44

        double precision, dimension(8,8), intent(OUT) :: k_el   !element stiffness matrix               
    
        double precision, dimension(4), intent(IN) :: x_e  !rectangular element (element coordinates)
        double precision, dimension(4), intent(IN) :: y_e  !rectangular element (element coordinates)

        double precision, dimension(4,4) :: k_u                         !see the shape2D_rectangle.f90 
        double precision, dimension(4,4) :: k_v                         !       -||- 
        double precision, dimension(4,4) :: k_uv                        !       -||- 
        double precision, dimension(4,4) :: k_vu                        !       -||- 
        
        double precision, dimension(4,4) :: FxFx                        !       -||- 
        double precision, dimension(4,4) :: FyFy                        !       -||- 
        double precision, dimension(4,4) :: FxFy                        !       -||- 
        double precision, dimension(4,4) :: FyFx                        !       -||-
        

        !##########-cubic elasticity-###################################################################
        !13.03.2017
                C_11 = C_11_i*(cos(theta_e)**4 + sin(theta_e)**4) + 2.0d+0*C_12_i*sin(theta_e)**2*cos(theta_e)**2 + C_44_i/4.0d+0*sin(2*theta_e)**2
                C_12 = C_12_i*(cos(theta_e)**4 + sin(theta_e)**4) + 2.0d+0*C_11_i*sin(theta_e)**2*cos(theta_e)**2 - C_44_i/4.0d+0*sin(2*theta_e)**2
                C_13 = 0.5d+0*C_11_i*sin(2*theta_e)*cos(2*theta_e) - 0.5d+0*C_12_i*sin(2*theta_e)*cos(2*theta_e) - 0.5d+0*C_44_i/4.0d+0*sin(4*theta_e)
                C_44 = 0.5d+0*C_11_i*sin(2*theta_e)**2 - 0.5d+0*C_12_i*sin(2*theta_e)**2 + C_44_i/4.0d+0*cos(2*theta_e)**2
                        C_44 = C_44*4.0d+0      !awkwardness due to classical definition of the shear components and tensorial shear strain
                                                !epsilon_xy = gamma_xy/2

        !##########-start stiffness matrix-##############################################################
        !this could be a subroutine, but can one subroutine have another subroutine inside
                                call shape2D_rectangle(FxFx, x_e,y_e,0)       !stiffness matrix
                                call shape2D_rectangle(FyFy, x_e,y_e,-1)      !stifness matrix
                                call shape2D_rectangle(FxFy, x_e,y_e,1)
                                FyFx = TRANSPOSE(FxFy)
                                
                                k_u = C_11*FxFx + C_44/4.0d+0*FyFy + C_13/1.0d+0*(FxFy+FyFx)
                                k_v = C_11*FyFy + C_44/4.0d+0*FxFx - C_13/1.0d+0*(FxFy+FyFx)
                                k_uv = C_12*FxFy + C_44/4.0d+0*FyFx + C_13/1.0d+0*(FxFx-FyFy)
                                k_vu = C_12*FyFx + C_44/4.0d+0*FxFy + C_13/1.0d+0*(FxFx-FyFy)
        
        !8x8 matrix 
!MATRIX ASSEMBLY
                DO,i=1,4,1      !4 nodes 
                DO,j=1,4,1      !4 nodes 
                        k_el(2*i-1,2*j-1) = k_u(i,j)    ![1,1]  
                        k_el(2*i-0,2*j-0) = k_v(i,j)    ![2,2]  
                        k_el(2*i-1,2*j-0) = k_uv(i,j)   ![1,2]  
                        k_el(2*i-0,2*j-1) = k_vu(i,j)   ![2,1]  
                END DO 
        END DO 
               
        !##########-end stiffness matrix-################################################################

END SUBROUTINE
