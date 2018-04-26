TFEM = [1,2,3,4,9,10,7,8
        3,4,5,6,11,12,9,10
        7,8,9,10,15,16,13,14
        9,10,11,12,17,18,15,16];
        
x_e = [-1 -1 0 0 0 0 -1 -1
       0 0 1 1 1 1 0 0 
       -1 -1 0 0 0 0 -1 -1
       0 0 1 1 1 1 0 0];
       
y_e = [-1 -1 -1 -1 0 0 0 0 
       -1 -1 -1 -1 0 0 0 0 
       0 0 0 0 1 1 1 1 
       0 0 0 0 1 1 1 1];
       
perm_element = [];
perm_dof2 = [];
indx_i = [];
indx_j = [];

for i=1:18 %nodes
[a_element,b_node] = find(TFEM==i);
A = sortrows([a_element,b_node]);
perm_element = [perm_element;A(:,1)];
perm_dof2 = [perm_dof2;A(:,2)];
indx_i = [indx_i;diag(x_e(A(:,1),A(:,2)))];
indx_j = [indx_j;diag(y_e(A(:,1),A(:,2)))];
end

%element --> node
%1 --> 5
%2 --> 7
%3 --> 3
%4 --> 1
				
perm_dof1(find(perm_element==1)) = 5;
perm_dof1(find(perm_element==2)) = 7;
perm_dof1(find(perm_element==3)) = 3;
perm_dof1(find(perm_element==4)) = 1;
perm_dof1 = perm_dof1';
