!!integral of shape functions 
!!Created by: Jakub Mikula
!!Created on: 31.10.2015
!!Updated on: 12.03.2016
!       contiguous arrays, all symmetric, : on the left, or rewritten

!##############################################################
!NOTES:
!
!##############################################################
!EXAMPLES:
!
!##############################################################
!shape function f = f(x,y)
!THE INTEGRALS ARE EXACTLY AS WRITTEN NEXT TO EACH OPTION, NO CHANGE IN SIGNS 

subroutine shape2D_rectangle(k_e, x_e,y_e,option)
implicit NONE

!declare variables
        double precision, dimension(4,4), intent(OUT) :: k_e
        double precision, dimension(4), intent(IN) :: x_e
        double precision, dimension(4), intent(IN) :: y_e
        
!        double precision, dimension(4) :: k_ex
!        double precision, dimension(4) :: k_ey
        
        integer :: option
        integer :: i
        integer :: j
        
        double precision :: al,bl
        
!options
        !option = 0     !\int \varphi_x \varphi_x dxdy                          !stiffness matrix 
        !option = -1 !\int \varphi_y \varphi_y dxdy                     !stiffness matrix 
        !option = 1 !\int \varphi_x \varphi_y dxdy      for \int \varphi_y \varphi_x dxdy TRANSPOSE the matrix 
        !option = 2 !\int \varphi \varphi dxdy                                  !mass matrix
        !option = 3 !\int \varphi_x \varphi dxdy

!rectangle dimensions 
        al = (x_e(2) - x_e(1))/2.0d+0   !why there is /2 ??? - because the length of the element is 2a and 2b, the area is 2ax2b
        if (al == 0) then
                al = (x_e(3) - x_e(2))/2.0d+0
        end if 
        
        bl = (y_e(2) - y_e(1))/2.0d+0 
        if (bl == 0) then 
                bl = (y_e(3) - y_e(2))/2.0d+0
    end if

    if (al < 1.0d-5 .OR. bl < 1.0d-5) THEN
        write(*,*) 'negative dimensions of the element <in shape2D_rectangle.f90>'
        stop
    end if


!element stiffness matrix
!shape_function_integration_bilinear.mw
!always access slices as v(:,:,1), ... etc., That way the stride is contiguous!
!always put the colon on the left
IF (option == 0) THEN               !\int \varphi_x \varphi_x dxdy                              !stiffness matrix 
        
                !k_e(1,:) = (/16,-16,-8,8/)
                !k_e(2,:) = (/-16,16,8,-8/)
                !k_e(3,:) = (/-8,8,16,-16/)
                !k_e(4,:) = (/8,-8,-16,16/)
                
        !contiguous (colon on the left!):
                k_e(:,1) = (/16,-16,-8,8/)
                k_e(:,2) = (/-16,16,8,-8/)
                k_e(:,3) = (/-8,8,16,-16/)
                k_e(:,4) = (/8,-8,-16,16/)

                k_e = bl/3.0/al*k_e/16.0 

ELSE IF (option == -1) THEN         !\int \varphi_y \varphi_y dxdy                      !stiffness matrix
        
                k_e(:,1) = (/16,8,-8,-16/)
                k_e(:,2) = (/8,16,-16,-8/)
                k_e(:,3) = (/-8,-16,16,8/)
                k_e(:,4) = (/-16,-8,8,16/)
                
                k_e = al/3.0/bl*k_e/16.0
        
ELSE IF (option == 1) THEN          !\int \varphi_x \varphi_y dxdy      for \int \varphi_y \varphi_x dxdy TRANSPOSE the matrix 

                !k_e(1,:) = (/4,4,-4,-4/)/16.0
                !k_e(2,:) = (/-4,-4,4,4/)/16.0 
                !k_e(3,:) = (/-4,-4,4,4/)/16.0 
                !k_e(4,:) = (/4,4,-4,-4/)/16.0 

        !contiguous (colon on the left!): 
                k_e(:,1) = (/4,-4,-4,4/)/16.0
                k_e(:,2) = (/4,-4,-4,4/)/16.0
                k_e(:,3) = (/-4,4,4,-4/)/16.0
                k_e(:,4) = (/-4,4,4,-4/)/16.0
        
ELSE IF (option == 2) THEN              !\int \varphi \varphi dxdy                                      !mass matrix

        !k_e(1,:) = (/4,2,1,2/)
        !k_e(2,:) = (/2,4,2,1/)
        !k_e(3,:) = (/1,2,4,2/)
        !k_e(4,:) = (/2,1,2,4/)
       
        !contiguous (colon on the left!):
        k_e(:,1) = (/4,2,1,2/)
        k_e(:,2) = (/2,4,2,1/)
        k_e(:,3) = (/1,2,4,2/)
        k_e(:,4) = (/2,1,2,4/)
 
        k_e = al*bl/9.0*k_e

ELSE IF (option == 3) THEN          !\int \varphi_x \varphi dxdy

        !k_e(1,:) = (/-2,2,1,-1/)
        !k_e(2,:) = (/-2,2,1,-1/)
        !k_e(3,:) = (/-1,1,2,-2/)
        !k_e(4,:) = (/-1,1,2,-2/)

        !contiguou (colon on the left!):
        k_e(:,1) = (/-2,-2,-1,-1/) 
        k_e(:,2) = (/2,2,1,1/)
        k_e(:,3) = (/1,1,2,2/)
        k_e(:,4) = (/-1,-1,-2,-2/)
        
        k_e = bl/6.0*k_e
        
ELSE IF (option == -3) THEN

        !k_e(1,:) = (/-2,-1,1,2/)
        !k_e(2,:) = (/-1,-2,2,1/)
        !k_e(3,:) = (/-1,-2,2,1/)
        !k_e(4,:) = (/-2,-1,1,2/)
        
        !contigous: 
        k_e(:,1) = (/-2,-1,-1,-2/)
        k_e(:,2) = (/-1,-2,-2,-1/)
        k_e(:,3) = (/1,2,2,1/)
        k_e(:,4) = (/2,1,1,2/)

        k_e = al/6.0*k_e
        
END IF

end subroutine
